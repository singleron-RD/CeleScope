
import numpy as np
import pysam
import sys
import subprocess

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT
from celescope.snp.__init__ import PANEL


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
        self.match_barcode_list, self.n_cell = utils.get_barcode_from_match_dir(args.match_dir)
        self.match_barcode = set(self.match_barcode_list)

        if args.panel:
            self.gene_list = utils.get_gene_region_from_bed(args.panel)[0]
            self.n_gene = len(self.gene_list)
        else:
            self.gene_list, self.n_gene = utils.read_one_col(args.gene_list)

        if not self.gene_list:
            sys.exit("You must provide either --panel or --gene_list!")

        self.count_dict = utils.genDict(dim=3, valType=int)

        self.add_metric(
            name="Number of Target Genes",
            value=self.n_gene,
        )
        self.add_metric(
            name="Number of Cells",
            value=self.n_cell,
        )

        # out file
        self.out_bam_file = f'{self.out_prefix}_filtered.bam'
        self.out_bam_file_sorted = f'{self.out_prefix}_filtered_sorted.bam'

    @utils.add_log
    def read_bam_write_filtered(self):
        sam_temp = f'{self.out_bam_file}.temp'
        with pysam.AlignmentFile(self.args.bam, "rb") as reader:
            header = reader.header.to_dict()
            # add RG to header
            if self.args.add_RG:
                header['RG'] = []
                for barcode in self.match_barcode_list:
                    header['RG'].append({
                        'ID': barcode,
                        'SM': barcode,
                    })
            with pysam.AlignmentFile(sam_temp, "w", header=header) as writer:
                for record in reader:
                    try:
                        gene_name = record.get_tag('GN')
                    except KeyError:
                        continue
                    # compatible with 10X bam
                    try:
                        barcode = record.get_tag('CB')
                        UMI = record.get_tag('UB')
                    except KeyError:
                        continue
                    if barcode in self.match_barcode and gene_name in self.gene_list:
                        if self.args.add_RG:
                            record.set_tag(tag='RG', value=record.get_tag('CB'), value_type='Z')
                        writer.write(record)
                    self.count_dict[barcode][gene_name][UMI] += 1

            cmd = f'samtools view -b {sam_temp} -o {self.out_bam_file}; rm {sam_temp}'
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def parse_count_dict_add_metrics(self):
        total_reads = 0
        enriched_reads = 0
        enriched_reads_in_cells = 0
        enriched_reads_per_cell_list = []

        for barcode in self.count_dict:
            cell_enriched_read = 0
            for gene_name in self.count_dict[barcode]:
                gene_read = sum(self.count_dict[barcode][gene_name].values())
                total_reads += gene_read
                if gene_name in self.gene_list:
                    enriched_reads += gene_read
                    if barcode in self.match_barcode:
                        enriched_reads_in_cells += gene_read
                        cell_enriched_read += gene_read
            if barcode in self.match_barcode:
                enriched_reads_per_cell_list.append(cell_enriched_read)

        self.parse_count_dict_add_metrics.logger.debug(
            f'enriched_reads_per_cell_list: '
            f'{sorted(enriched_reads_per_cell_list)}'
            f'len: {len(enriched_reads_per_cell_list)}'
        )

        valid_enriched_reads_per_cell_list = [cell for cell in enriched_reads_per_cell_list if cell > 0]
        n_valid_cell = len(valid_enriched_reads_per_cell_list)
        self.add_metric(
            name="Number of Valid Cells",
            value=n_valid_cell,
            total=self.n_cell,
        )
        self.add_metric(
            name="Enriched Reads",
            value=enriched_reads,
            total=total_reads,
        )
        self.add_metric(
            name="Enriched Reads in Cells",
            value=enriched_reads_in_cells,
            total=total_reads,
        )
        self.add_metric(
            name="Median Enriched Reads per Valid Cell",
            value=np.median(valid_enriched_reads_per_cell_list),
        )

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
    with Target_metrics(args, display_title='Target Enrichment') as runner:
        runner.run()


def get_opts_target_metrics(parser, sub_program):
    parser.add_argument("--gene_list", help=HELP_DICT['gene_list'])
    parser.add_argument("--panel", help=HELP_DICT['panel'], choices=list(PANEL))
    if sub_program:
        parser.add_argument("--bam", help='Input bam file', required=True)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument(
            '--add_RG', help='Add tag read group: RG. RG is the same as CB(cell barcode)', action='store_true')
        parser = s_common(parser)
