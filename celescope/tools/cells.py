from collections import defaultdict
import shutil
import os

import numpy as np
import pandas as pd

from celescope.tools.__init__ import FILTERED_MATRIX_DIR_SUFFIX, RAW_MATRIX_DIR_SUFFIX, COUNTS_FILE_NAME 
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.tools.matrix import CountMatrix
from celescope.tools.emptydrop_cr import get_plot_elements

class Cells_metrics(Step):
    @utils.add_log
    def add_cells_metrics(self, n_cells, fraction_reads_in_cells, mean_used_reads_per_cell, median_umi_per_cell,
        total_genes, median_genes_per_cell, saturation):
        self.add_metric(
            name='Estimated Number of Cells',
            value=n_cells,
            help_info='the number of barcodes considered as cell-associated.'
        )

        self.add_metric(
            name='Fraction Reads in Cells',
            value=fraction_reads_in_cells,
            value_type='fraction',
            help_info='the fraction of uniquely-mapped-to-transcriptome reads with cell-associated barcodes'
        )

        self.add_metric(
            name='Mean Used Reads per Cell',
            value=mean_used_reads_per_cell,
            help_info='the number of uniquely-mapped-to-transcriptome reads per cell-associated barcode'
        )

        self.add_metric(
            name='Median UMI per Cell',
            value=median_umi_per_cell,
            help_info='the median number of UMI counts per cell-associated barcode'
        )

        self.add_metric(
            name='Total Genes',
            value=total_genes,
            help_info='the number of genes with at least one UMI count in any cell'
        )

        self.add_metric(
            name='Median Genes per Cell',
            value=median_genes_per_cell,
            help_info='the median number of genes detected per cell-associated barcode'
        )

        self.add_metric(
            name='Saturation',
            value=saturation,
            value_type='fraction',
            help_info='the fraction of read originating from an already-observed UMI. '
        )

    def run(self):
        pass

class Cells(Cells_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        self.raw_matrix = f'{self.outs_dir}/{RAW_MATRIX_DIR_SUFFIX}'
        self.old_filtered_matrix = f'{self.outs_dir}/{FILTERED_MATRIX_DIR_SUFFIX}'
        self.counts_file = f'{self.outs_dir}/{COUNTS_FILE_NAME}'

        # out
        self.filter_matrix = f'{self.outdir}/{FILTERED_MATRIX_DIR_SUFFIX}'
        self.default_filter_matrix = f'{self.outs_dir}/default_{FILTERED_MATRIX_DIR_SUFFIX}'

        self.outs = [self.filter_matrix]

    @utils.add_log
    def force_cells(self):
        if self.args.force_cells <= 0:
            return None, []
        raw = CountMatrix.from_matrix_dir(self.raw_matrix)
        df_counts = pd.read_csv(self.counts_file, index_col=0, header=0, sep='\t')
        bcs = list(df_counts.head(self.args.force_cells).index)
        filtered = raw.slice_matrix_bc(bcs)
        if not os.path.exists(self.default_filter_matrix):
            shutil.move(self.old_filtered_matrix, self.default_filter_matrix)
        filtered.to_matrix_dir(self.filter_matrix)
        return filtered, bcs

    @utils.add_log
    def metrics_report(self, filtered, bcs):

        n_cells = len(bcs)

        df_counts = pd.read_csv(self.counts_file, index_col=0, header=0, sep='\t')
        reads_total = df_counts['countedU'].sum()
        reads_cell = df_counts.loc[bcs, 'countedU'].sum()
        fraction_reads_in_cells = float(reads_cell / reads_total)
        mean_used_reads_per_cell = int(reads_cell // len(bcs))
        median_umi_per_cell = int(df_counts.loc[bcs, 'UMI'].median())

        matrix = filtered.get_matrix()
        gene_index, bc_index = matrix.nonzero()
        total_genes = len(set(gene_index))
        bc_gene = defaultdict(set)
        for gene, bc in zip(gene_index, bc_index):
            bc_gene[bc].add(gene)
        gene_per_cell = [len(v) for v in bc_gene.values()]            
        median_genes_per_cell = int(np.median(gene_per_cell))

        saturation = self.old_step_dict['metrics']['Saturation'] / 100
        df_counts.loc[:, 'mark'] = 'UB'
        df_counts.loc[bcs, 'mark'] = 'CB'
        df_counts.to_csv(self.counts_file, sep='\t', index=True)
        self.add_cells_metrics(n_cells, fraction_reads_in_cells, mean_used_reads_per_cell, median_umi_per_cell, total_genes, median_genes_per_cell, saturation)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.counts_file))

    def run(self):
        bcs = []
        if self.args.force_cells > 0:
            filtered, bcs = self.force_cells()
        if bcs:
            self.metrics_report(filtered, bcs)

def cells(args):
    with Cells(args) as runner:
        runner.run()


def get_opts_cells(parser, sub_program=True):
    if sub_program:
        parser.add_argument(
            '--force_cells',
            help='Force to use this number of cells.',
            default=0,
            type=int,
        )
        parser.add_argument(
            '--root_dir',
            help='Root directory of CeleScope',
            default='./',
        )
        parser = s_common(parser)