"""
count step
"""

import os
import random
import sys

import numpy as np
import pandas as pd

from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.tools.__init__ import (FILTERED_MATRIX_DIR_SUFFIX, RAW_MATRIX_DIR_SUFFIX)
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools.emptydrop_cr.cell_calling_3 import cell_calling_3
from celescope.tools.step import Step, s_common
from celescope.rna.mkref import Mkref_rna
from celescope.tools.plotly_plot import Line_plot
from celescope.tools.matrix import CountMatrix, ROW, COLUMN
from celescope.tools import reference

TOOLS_DIR = os.path.dirname(__file__)
random.seed(0)
np.random.seed(0)

# downsample.csv
READ_FRACTION = 'read_fraction'
MEDIAN_GENE_NUMBER = 'median_gene_number'
SATURATION = 'saturation'

# Plot axis title in HTML
X_TITLE = 'Read Fraction'
SATURATION_Y_TITLE = 'Sequencing Saturation(%)'
MEDIAN_GENE_Y_TITLE = 'Median Genes per Cell'


class Count(Step):
    """
    ## Features
    - Cell-calling: Distinguish cell barcodes from background barcodes. 
    - Generate expression matrix.
    ## Output
    - `{sample}_raw_feature_bc_matrix` The expression matrix of all detected barcodes in [Matrix Market Exchange Formats](
        https://math.nist.gov/MatrixMarket/formats.html). 
    - `{sample}_filtered_feature_bc_matrix` The expression matrix of cell barcodes in Matrix Market Exchange Formats. 
    - `{sample}_count_detail.txt.gz` 4 columns: 
        - barcode  
        - gene ID  
        - UMI count  
        - read_count  
    - `{sample}_counts.txt` 6 columns:
        - Barcode: barcode sequence
        - read: read count of each barcode
        - UMI: UMI count for each barcode
        - geneID: gene count for each barcode
        - mark: cell barcode or backgound barcode.
            `CB` cell  
            `UB` background  
    - `{sample}_downsample.tsv` Subset a fraction of reads and calculate median gene number and sequencing saturation.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.force_cell_num = args.force_cell_num
        self.cell_calling_method = args.cell_calling_method
        self.expected_cell_num = int(args.expected_cell_num)

        # set
        gtf_file = Mkref_rna.parse_genomeDir(args.genomeDir)['gtf']
        gp = reference.GtfParser(gtf_file)
        gp.get_id_name()
        self.features = gp.get_features()
        self.downsample_dict = {}

        # output files
        self.marked_count_file = f'{self.outdir}/{self.sample}_counts.txt'
        self.raw_matrix_dir = f'{self.outdir}/{self.sample}_{RAW_MATRIX_DIR_SUFFIX[0]}'
        self.cell_matrix_dir = f'{self.outdir}/{self.sample}_{FILTERED_MATRIX_DIR_SUFFIX[0]}'
        self.downsample_file = f'{self.outdir}/{self.sample}_downsample.tsv'

    def get_df_line(self):
        df_line = pd.read_csv(self.downsample_file, sep="\t", header=0)
        df_line.rename(columns={
            READ_FRACTION: X_TITLE,
            MEDIAN_GENE_NUMBER: MEDIAN_GENE_Y_TITLE,
            SATURATION: SATURATION_Y_TITLE,
        }, inplace=True)

        return df_line

    @staticmethod
    def get_df_cell(df, cell_bc):
        """
        reset index in a must. otherwise the index remains the same with the raw matrix
        """
        df_cell = df.loc[df.index.isin(cell_bc, level=0)]
        df_cell = df_cell.reset_index().set_index([COLUMN, ROW])
        return df_cell

    @utils.add_log
    def run(self):
        # index use less than >50% mem than int index
        df = pd.read_table(self.args.count_detail, header=0, dtype={'geneID':str}, index_col=[0,1])
        total_reads = sum(df['read'])
        self.write_sparse_matrix(df, self.raw_matrix_dir)
        df_bc = self.get_df_bc(df)
        # call cells
        cell_bc, _threshold = self.cell_calling(df_bc)
        # write marked_df_sum
        self.write_marked_df_bc(df_bc, cell_bc)
        # export cell matrix
        df_cell = self.get_df_cell(df, cell_bc)
        del df
        self.write_sparse_matrix(df_cell, self.cell_matrix_dir)
        total_cell_gene = len(df_cell.index.levels[1])
        df_bc_cell = self.get_df_bc(df_cell)
        del df_bc
        # downsample
        downsample_dict = self.downsample(df_cell)
        # metrics
        self.add_count_metrics(total_reads, total_cell_gene, cell_bc, df_bc_cell, downsample_dict)
        # plot
        self.add_plot_data()

    def add_count_metrics(self, total_reads, total_cell_gene, cell_bc, df_sum_cell, downsample_dict):
        n_cells = len(cell_bc)
        self.add_metric(
            name='Estimated Number of Cells',
            value=n_cells,
            help_info=f'the number of barcodes considered as cell-associated. The cell calling method used for this sample is {self.cell_calling_method}.'
        )

        cell_reads = sum(df_sum_cell['read'])
        fraction_reads_in_cells = round(float(cell_reads) / total_reads * 100, 2)
        self.add_metric(
            name='Fraction Reads in Cells',
            value=fraction_reads_in_cells,
            display=f'{fraction_reads_in_cells}%',
            help_info='the fraction of uniquely-mapped-to-transcriptome reads with cell-associated barcodes'
        )

        try:
            valid_read_number = self.get_slot_key(
                slot='metrics',
                step_name='barcode',
                key='Valid Reads',
            )
        except KeyError:
            self.add_count_metrics.logger.warning('Will not output `Mean Reads per Cell`')
        else:
            mean_reads_per_cell = int(valid_read_number / n_cells)
            self.add_metric(
                name='Mean Reads per Cell',
                value=mean_reads_per_cell,
                help_info='the number of valid reads divided by the estimated number of cells'
            )

        median_umi_per_cell = int(np.median(df_sum_cell['UMI']))
        self.add_metric(
            name='Median UMI per Cell',
            value=median_umi_per_cell,
            help_info='the median number of UMI counts per cell-associated barcode'
        )

        self.add_metric(
            name='Total Genes',
            value=total_cell_gene,
            help_info='the number of genes with at least one UMI count in any cell'
        )

        self.add_metric(
            name='Median Genes per Cell',
            value=int(downsample_dict[MEDIAN_GENE_NUMBER][-1]),
            help_info='the median number of genes detected per cell-associated barcode'
        )

        saturation = downsample_dict[SATURATION][-1]
        self.add_metric(
            name='Saturation',
            value=saturation,
            display=f'{saturation}%',
            help_info='the fraction of UMI originating from an already-observed UMI. '
        )


    def add_plot_data(self):

        df_line = self.get_df_line()

        line_saturation = Line_plot(df_line=df_line,title="Sequencing Saturation",x_title=X_TITLE,
                                    y_title=SATURATION_Y_TITLE, y_range=[0,100], section=False).get_plotly_div()
        self.add_data(line_saturation=line_saturation)
        line_median = Line_plot(df_line=df_line,title="Median Genes per Cell", x_title=X_TITLE,
                                y_title=MEDIAN_GENE_Y_TITLE).get_plotly_div()
        self.add_data(line_median=line_median)

        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.marked_count_file))


    
    @utils.add_log
    def cell_calling(self, df_sum):
        """
        Returns:
            cell_bc: set
            UMI_threshold: int
        """
        cell_calling_method = self.cell_calling_method

        if (self.force_cell_num is not None) and (self.force_cell_num != 'None'):
            cell_bc, UMI_threshold = self.force_cell(df_sum)
        elif cell_calling_method == 'auto':
            cell_bc, UMI_threshold = self.auto_cell(df_sum)
        elif cell_calling_method == 'EmptyDrops_CR':
            cell_bc, UMI_threshold = self.emptydrop_cr_cell(df_sum)
        cell_bc = set(cell_bc)
        return cell_bc, UMI_threshold

    @utils.add_log
    def force_cell(self, df_sum):
        force_cell_num = int(self.force_cell_num)

        df_barcode_count = df_sum.groupby(
            ['UMI']).size().reset_index(
            name='barcode_counts')
        sorted_df = df_barcode_count.sort_values("UMI", ascending=False)
        sorted_df["barcode_cumsum"] = sorted_df["barcode_counts"].cumsum()
        num_points = sorted_df.shape[0]
        for i in range(num_points):
            if sorted_df.iloc[i, :]["barcode_cumsum"] >= force_cell_num:
                index_low = i - 1
                index_high = i
                break
        df_sub = sorted_df.iloc[index_low: index_high + 1, :]
        distance = abs(df_sub["barcode_cumsum"] - force_cell_num)
        actual_index = np.argmin(distance)
        threshold = df_sub.iloc[actual_index, :]['UMI']
        cell_bc = Count.get_cell_bc(df_sum, threshold, col='UMI')
        return cell_bc, threshold

    @staticmethod
    def find_threshold(df_sum, idx):
        return int(df_sum.iloc[idx - 1, df_sum.columns == 'UMI'])

    @staticmethod
    def get_cell_bc(df_sum, threshold, col='UMI'):
        return list(df_sum[df_sum[col] >= threshold].index)

    @utils.add_log
    def auto_cell(self, df_sum):
        idx = int(self.expected_cell_num * 0.01)
        barcode_number = df_sum.shape[0]
        idx = int(min(barcode_number, idx))
        if idx == 0:
            sys.exit("cell number equals zero!")
        # calculate read counts threshold
        threshold = int(Count.find_threshold(df_sum, idx) * 0.1)
        threshold = max(1, threshold)
        cell_bc = Count.get_cell_bc(df_sum, threshold)

        return cell_bc, threshold

    @utils.add_log
    def emptydrop_cr_cell(self, df_sum):
        cell_bc, initial_cell_num = cell_calling_3(self.raw_matrix_dir, self.expected_cell_num)
        threshold = Count.find_threshold(df_sum, initial_cell_num)
        return cell_bc, threshold

    @staticmethod
    def get_df_bc(df, sort_by='UMI'):
        '''
        Returns:
            df with 3 cols sort by UMI(default)
                read: int
                UMI: int
                geneID: str
        '''

        df_bc = df.groupby('Barcode').agg({
            'UMI': 'sum',
            'read': 'sum',
        })
        df_bc[ROW] = df.groupby(level=[0]).size()
        return df_bc.sort_values(sort_by, ascending=False)

    def write_marked_df_bc(self, df_sum, cell_bc):
        df_sum.loc[:, 'mark'] = 'UB'
        df_sum.loc[df_sum.index.isin(cell_bc), 'mark'] = 'CB'
        df_sum.to_csv(self.marked_count_file, sep='\t')

    @utils.add_log
    def write_sparse_matrix(self, df, matrix_dir):

        count_matrix = CountMatrix.from_dataframe(df, self.features, value="UMI")
        count_matrix.to_matrix_dir(matrix_dir)

    @utils.add_log
    def cell_summary(self, df, cell_bc):

        df.loc[:, 'mark'] = 'UB'
        df.loc[df['Barcode'].isin(cell_bc), 'mark'] = 'CB'
        CB_total_Genes = df.loc[df['mark'] == 'CB', 'geneID'].nunique()
        CB_reads_count = df.loc[df['mark'] == 'CB', 'count'].sum()
        reads_mapped_to_transcriptome = df['count'].sum()
        return(CB_total_Genes, CB_reads_count, reads_mapped_to_transcriptome)

    @staticmethod
    def sub_sample(fraction, df_cell, new_df, cell_read_index):
        """
        saturation = 1 - n_deduped_reads / n_reads

        n_deduped_reads = Number of unique (valid cell-barcode, valid UMI, gene) combinations among confidently mapped reads. If a UMI has two reads mapped to one gene, but mapped to different locations of the gene, these two reads are
        considered unique.
        n_reads = Total number of (confidently mapped, valid cell-barcode, valid UMI) reads.
        Args:
            fration: subsmaple fration
            df_cell: in cell df with (Barcode geneID UMI count) 
            new_df: two columns: [original df_cell index, duplicate count]
            cell_read_index: df_cell repeat index
        """

        cell_read = new_df['dcount'].sum()
        frac_n_read = int(cell_read * fraction)
        subsample_read_index = cell_read_index[:frac_n_read]
        index_dedup, counts = np.unique(subsample_read_index, return_counts=True)
        n_count_once = np.sum(counts == 1)
        read_total = frac_n_read
        saturation = round((1 - n_count_once / read_total) * 100, 2)
        df_cell_subsample = df_cell.loc[new_df.loc[index_dedup,]['oindex']]
        geneNum_median = float(df_cell_subsample.groupby(
        'Barcode').agg({'geneID': 'nunique'}).median())

        return saturation, geneNum_median

    @utils.add_log
    def downsample(self, df_cell):
        """saturation and median gene
        Returns:
            downsample dict             
                READ_FRACTION: float
                SATURATION: float
                MEDIAN_GENE_NUMBER: float
        """
        oindexList = []
        dcountList = []
        df_cell = df_cell.reset_index()

        for row in df_cell.itertuples():
            duplicates = row.duplicate.split(',')
            ds = [int(d) for d in duplicates]
            for d in ds:
                oindexList.append(row.Index)
                dcountList.append(d)

        new_df = pd.DataFrame({'oindex':oindexList, 'dcount':dcountList})

        cell_read_index = np.array(new_df.index.repeat(new_df['dcount']), dtype='int32')
        np.random.shuffle(cell_read_index)

        downsample_dict = {
            READ_FRACTION: [0],
            SATURATION: [0],
            MEDIAN_GENE_NUMBER: [0],
        }

        for fraction in np.arange(0.1, 1.1, 0.1):
            saturation, geneNum_median = Count.sub_sample(fraction, df_cell, new_df, cell_read_index)
            fraction = round(fraction,1)
            downsample_dict[READ_FRACTION].append(fraction)
            downsample_dict[SATURATION].append(saturation)
            downsample_dict[MEDIAN_GENE_NUMBER].append(geneNum_median)

        df_downsample = pd.DataFrame(downsample_dict, columns=[READ_FRACTION, MEDIAN_GENE_NUMBER, SATURATION,])
        df_downsample.to_csv(self.downsample_file, index=False, sep='\t')
        return downsample_dict


@utils.add_log
def count(args):
    with Count(args, display_title="Cells") as runner:
        runner.run()


def get_opts_count(parser, sub_program):
    parser.add_argument('--genomeDir', help='Required. Genome directory.')
    parser.add_argument('--expected_cell_num', help='Default `3000`. Expected cell number.', default=3000)
    parser.add_argument(
        '--cell_calling_method',
        help=HELP_DICT['cell_calling_method'],
        choices=['auto', 'EmptyDrops_CR'],
        default='EmptyDrops_CR',
    )
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--count_detail', help='Required. File from featureCounts.', required=True)
        parser.add_argument(
            '--force_cell_num',
            help='Default `None`. Force the cell number to be this number. ',
        )
