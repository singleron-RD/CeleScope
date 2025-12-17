import numpy as np
import pandas as pd
import pysam
import sys
import subprocess
from collections import defaultdict

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT
from celescope.snp.__init__ import PANEL


def get_genes(args) -> set:
    if args.bed:
        genes, _ = utils.get_gene_region_from_bed(args.bed)
    elif args.panel:
        bed = utils.get_bed_file_path(args.panel)
        genes, _ = utils.get_gene_region_from_bed(bed)
    elif args.gene_list:
        genes, _ = utils.read_one_col(args.gene_list)
    else:
        sys.exit("You must provide one of --panel, --bed, --gene_list!")

    return set(genes)


class Target_metrics(Step):
    """
    ## Features
    - Filter bam file
        - Filter reads that are not cell-associated.
        - Filter reads that are not mapped to target genes.

    - Collect enrichment metrics.

    ## Output
    - `filtered.bam` BAM file after filtering. Reads that are not cell-associated or not mapped to target genes are filtered.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # set
        self.match_barcode_list, self.n_cell = utils.get_barcode_from_match_dir(
            args.match_dir
        )
        self.match_barcode = set(self.match_barcode_list)

        self.genes = get_genes(args)

        self.count_dict = utils.nested_defaultdict(dim=3, valType=int)
        self.dup_dict = utils.nested_defaultdict(dim=4, valType=int)
        self.used_dict = defaultdict(int)

        self.add_metric(
            name="Number of Target Genes",
            value=len(self.genes),
        )
        self.add_metric(
            name="Number of Cells",
            value=self.n_cell,
        )

        # out file
        self.out_bam_file = f"{self.out_prefix}_filtered.bam"
        self.out_bam_file_sorted = f"{self.out_prefix}_filtered_sorted.bam"
        self.gene_count_tsv = f"{self.out_prefix}_gene_count.tsv"

    @utils.add_log
    def read_bam_write_filtered(self):
        """
        for each (barcode,UMI,reference_name,reference_start), keep at most max_duplicate_reads
        """
        sam_temp = f"{self.out_bam_file}.temp"
        max_duplicate = self.args.max_duplicate
        with pysam.AlignmentFile(self.args.bam, "rb") as reader:
            header = reader.header.to_dict()
            # add RG to header
            if self.args.add_RG:
                header["RG"] = []
                for barcode in self.match_barcode_list:
                    header["RG"].append(
                        {
                            "ID": barcode,
                            "SM": barcode,
                        }
                    )
            with pysam.AlignmentFile(sam_temp, "w", header=header) as writer:
                for record in reader:
                    try:
                        gene_name = record.get_tag("GN")
                    except KeyError:
                        continue
                    # compatible with tag bam
                    try:
                        barcode = record.get_tag("CB")
                        UMI = record.get_tag("UB")
                    except KeyError:
                        continue
                    self.count_dict[barcode][gene_name][UMI] += 1
                    if barcode in self.match_barcode and gene_name in self.genes:
                        rn, rs = record.reference_name, record.reference_start
                        self.dup_dict[barcode][UMI][rn][rs] += 1
                        if self.dup_dict[barcode][UMI][rn][rs] > max_duplicate:
                            continue
                        self.used_dict[barcode] += 1
                        if self.args.add_RG:
                            record.set_tag(
                                tag="RG", value=record.get_tag("CB"), value_type="Z"
                            )
                        writer.write(record)

            cmd = f"samtools view -b {sam_temp} -o {self.out_bam_file}; rm {sam_temp}"
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def parse_count_dict_add_metrics(self):
        total_reads = 0
        enriched_reads = 0
        enriched_reads_in_cells = 0
        enriched_reads_per_cell_list = []
        gene_count_dict = defaultdict(int)

        for barcode in self.count_dict:
            cell_enriched_read = 0
            for gene_name in self.count_dict[barcode]:
                gene_read = sum(self.count_dict[barcode][gene_name].values())
                total_reads += gene_read
                gene_count_dict[gene_name] += gene_read
                if gene_name in self.genes:
                    enriched_reads += gene_read
                    if barcode in self.match_barcode:
                        enriched_reads_in_cells += gene_read
                        cell_enriched_read += gene_read
            if barcode in self.match_barcode:
                enriched_reads_per_cell_list.append(cell_enriched_read)

        self.parse_count_dict_add_metrics.logger.debug(
            f"enriched_reads_per_cell_list: "
            f"{sorted(enriched_reads_per_cell_list)}"
            f"len: {len(enriched_reads_per_cell_list)}"
        )

        valid_enriched_reads_per_cell_list = [
            cell for cell in enriched_reads_per_cell_list if cell > 0
        ]
        n_valid_cell = len(valid_enriched_reads_per_cell_list)
        self.add_metric(
            name="Number of Valid Cells",
            value=n_valid_cell,
            total=self.n_cell,
            help_info="Cells with at least one read mapped to the target genes.",
        )
        self.add_metric(
            name="Enriched Reads",
            value=enriched_reads,
            total=total_reads,
            help_info="Reads mapped to the target genes.",
        )
        self.add_metric(
            name="Enriched Reads in Cells",
            value=enriched_reads_in_cells,
            total=total_reads,
            help_info="Reads mapped to the target genes and belong to cells(from scRNA match_dir).",
        )
        self.add_metric(
            name="Median Enriched Reads per Valid Cell",
            value=np.median(valid_enriched_reads_per_cell_list),
            help_info="Median number of reads mapped to the target genes per cell.",
        )

        self.add_metric(
            name="Median Used Reads per Valid Cell",
            value=np.median(list(self.used_dict.values())),
            help_info="Median number of reads per cell after deduplication. Deduplication is performed using the combination (CB, UB, reference_name, reference_start). For each unique combination, at most max_duplicate reads are retained.",
        )

        # write gene count tsv
        df = pd.DataFrame.from_dict(gene_count_dict, orient="index", columns=["count"])
        df.index.name = "gene_name"
        df = df.sort_values(by="count", ascending=False)
        df["pencent"] = round(df["count"] / total_reads * 100, 2)
        df.to_csv(self.gene_count_tsv, sep="\t")

    def run(self):
        self.read_bam_write_filtered()
        self.parse_count_dict_add_metrics()
        samtools_runner = utils.Samtools(
            self.out_bam_file,
            self.out_bam_file_sorted,
            self.args.thread,
            debug=self.debug,
        )
        samtools_runner.sort_bam()
        samtools_runner.index_bam()


@utils.add_log
def target_metrics(args):
    with Target_metrics(args, display_title="Target Enrichment") as runner:
        runner.run()


def get_opts_target_metrics(parser, sub_program):
    parser.add_argument("--gene_list", help=HELP_DICT["gene_list"])
    parser.add_argument("--panel", help=HELP_DICT["panel"], choices=list(PANEL))
    parser.add_argument("--bed", help="custom bed file.")
    parser.add_argument("--max_duplicate", type=int, default=5)
    if sub_program:
        parser.add_argument("--bam", help="Input bam file", required=True)
        parser.add_argument("--match_dir", help=HELP_DICT["match_dir"], required=True)
        parser.add_argument(
            "--add_RG",
            help="Add tag read group: RG. RG is the same as CB(cell barcode)",
            action="store_true",
        )
        parser = s_common(parser)
