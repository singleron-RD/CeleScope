
import numpy as np
import pysam
import sys

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
        self.match_barcode, self.n_cell = utils.read_barcode_file(args.match_dir)
        self.match_barcode = set(self.match_barcode)
        
        if (self.assay == "snp" and args.panel != '') :
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

        if self.debug:
            self.parse_count_dict_add_metrics.logger.debug(
                f'enriched_reads_per_cell_list: '
                f'{sorted(enriched_reads_per_cell_list)}'
                f'len: {len(enriched_reads_per_cell_list)}'
            )
        
        n_valid_cell = len([cell for cell in enriched_reads_per_cell_list if cell > 0])
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
            value=np.median(enriched_reads_per_cell_list),
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
    parser.add_argument("--gene_list", help=HELP_DICT['gene_list'])
    parser.add_argument("--panel",help = "The prefix of bed file, such as `lung_1`.",default = '')
    if sub_program:
        parser.add_argument("--bam", help='Input bam file', required=True)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser = s_common(parser)
