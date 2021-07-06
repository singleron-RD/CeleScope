
import numpy as np
import pysam

import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT


class Target_metrics(Step):
    """
    Features
    - Filter bam file
        - Filter reads that are not cell-associated.
        - Filter reads that are not mapped to target genes. 

    - Collect enrichment metrics.

    Output
    - `filtered.bam` BAM file after filtering.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # set
        self.match_barcode, _num = utils.read_barcode_file(args.match_dir)
        self.gene_list, self.n_gene = utils.read_one_col(args.gene_list)
        self.count_dict = utils.genDict(dim=3, valType=int)

        self.add_metric(
            name="Number of Target Genes",
            value=self.n_gene,
        )

        # out file
        self.out_bam_file = f'{self.out_prefix}_filtered.bam'
        self.out_bam_file_sorted = f'{self.out_prefix}_filtered_sorted.bam'

    @utils.add_log
    def read_bam_write_filtered(self):
        with pysam.AlignmentFile(self.args.bam, "rb") as reader:
            with pysam.AlignmentFile(self.out_bam_file, "wb", header=reader.header) as writer:
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
                        writer.write(record)
                    self.count_dict[barcode][gene_name][UMI] += 1

    @utils.add_log
    def parse_count_dict_add_metrics(self):
        total_UMIs = 0
        enriched_UMIs = 0
        enriched_UMIs_in_cells = 0
        enriched_UMIs_per_cell_list = []

        for barcode in self.count_dict:
            cell_enriched_UMI = 0
            for gene_name in self.count_dict[barcode]:
                gene_UMI = len(self.count_dict[barcode][gene_name])
                total_UMIs += gene_UMI
                if gene_name in self.gene_list:
                    enriched_UMIs += gene_UMI
                    if barcode in self.match_barcode:
                        enriched_UMIs_in_cells += gene_UMI
                        cell_enriched_UMI += gene_UMI

            if barcode in self.match_barcode:
                enriched_UMIs_per_cell_list.append(cell_enriched_UMI)

        self.add_metric(
            name="Total UMIs",
            value=total_UMIs,
        )

        self.add_metric(
            name="Enriched UMIs",
            value=enriched_UMIs,
            total=total_UMIs,
        )
        self.add_metric(
            name="Enriched UMIs in Cells",
            value=enriched_UMIs_in_cells,
            total=total_UMIs,
        )
        self.add_metric(
            name="Median Enriched UMIs per Cell",
            value=np.median(enriched_UMIs_per_cell_list),
        )

    def run(self):
        self.read_bam_write_filtered()
        self.parse_count_dict_add_metrics()
        utils.sort_bam(
            self.out_bam_file,
            self.out_bam_file_sorted,
            threads=self.thread,
        )
        utils.index_bam(self.out_bam_file_sorted)
        self.clean_up()


@utils.add_log
def target_metrics(args):
    step_name = "target_metrics"
    runner = Target_metrics(args, step_name)
    runner.run()


def get_opts_target_metrics(parser, sub_program):
    parser.add_argument("--gene_list", help=HELP_DICT['gene_list'], required=True)
    if sub_program:
        parser.add_argument("--bam", help='Input bam file', required=True)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser = s_common(parser)
