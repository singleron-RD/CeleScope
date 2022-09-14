"""
count step
"""

import os
import random
import sys
import unittest
from collections import defaultdict
from itertools import groupby

import numpy as np
import pandas as pd
import pysam

from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.tools.__init__ import (FILTERED_MATRIX_DIR_SUFFIX, RAW_MATRIX_DIR_SUFFIX)
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools.emptydrop_cr.cell_calling_3 import cell_calling_3
from celescope.tools.step import Step, s_common
from celescope.rna.mkref import Mkref_rna
from celescope.tools.plotly_plot import Line_plot
from celescope.tools.matrix import CountMatrix
from celescope.tools import reference

TOOLS_DIR = os.path.dirname(__file__)
random.seed(0)
np.random.seed(0)

# downsample.csv
READ_FRACTION = 'read_fraction'
MEDIAN_GENE_NUMBER = 'median_gene_number'
READ_SATURATION = 'read_saturation'
UMI_SATURATION = 'umi_saturation'

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
        - readcount: read count of each barcode
        - UMI2: read count with reads per UMI >= 2 for each barcode
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
        self.bam = args.bam

        # set
        gtf_file = Mkref_rna.parse_genomeDir(args.genomeDir)['gtf']
        gp = reference.GtfParser(gtf_file)
        gp.get_id_name()
        self.features = gp.get_features()
        self.downsample_dict = {}

        # output files
        self.count_detail_file = f'{self.outdir}/{self.sample}_count_detail.txt'
        self.marked_count_file = f'{self.outdir}/{self.sample}_counts.txt'
        self.raw_matrix_dir = f'{self.outdir}/{self.sample}_{RAW_MATRIX_DIR_SUFFIX[0]}'
        self.cell_matrix_dir = f'{self.outdir}/{self.sample}_{FILTERED_MATRIX_DIR_SUFFIX[0]}'
        self.downsample_file = f'{self.outdir}/{self.sample}_downsample.tsv'

    def get_df_line(self):
        df_line = pd.read_csv(self.downsample_file, sep="\t", header=0)
        df_line.rename(columns={
            READ_FRACTION: X_TITLE,
            MEDIAN_GENE_NUMBER: MEDIAN_GENE_Y_TITLE,
            UMI_SATURATION: SATURATION_Y_TITLE,
        }, inplace=True)

        return df_line

    def run(self):
        self.bam2table()
        df = pd.read_table(self.count_detail_file, header=0)

        # df_sum
        df_sum = Count.get_df_sum(df)

        # export all matrix
        self.write_sparse_matrix(df, self.raw_matrix_dir)

        # call cells
        cell_bc, _threshold = self.cell_calling(df_sum)

        # get cell stats
        CB_describe = self.get_cell_stats(df_sum, cell_bc)

        # export cell matrix
        df_cell = df.loc[df['Barcode'].isin(cell_bc), :]
        self.write_sparse_matrix(df_cell, self.cell_matrix_dir)
        (CB_total_Genes, CB_reads_count, reads_mapped_to_transcriptome) = self.cell_summary(
            df, cell_bc)

        # downsampling
        cell_bc = set(cell_bc)
        self.downsample(df_cell)

        # summary
        self.get_summary(CB_describe, CB_total_Genes,
                         CB_reads_count, reads_mapped_to_transcriptome)

        df_line = self.get_df_line()

        line_saturation = Line_plot(df_line=df_line,title="Sequencing Saturation",x_title=X_TITLE,
                                    y_title=SATURATION_Y_TITLE, y_range=[0,100], section=False).get_plotly_div()
        self.add_data(line_saturation=line_saturation)
        line_median = Line_plot(df_line=df_line,title="Median Genes per Cell", x_title=X_TITLE,
                                y_title=MEDIAN_GENE_Y_TITLE).get_plotly_div()
        self.add_data(line_median=line_median)

        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.marked_count_file))

    @staticmethod
    def correct_umi(umi_dict, percent=0.1):
        """
        Correct umi_dict in place.
        Args:
            umi_dict: {umi_seq: umi_count}
            percent: if hamming_distance(low_seq, high_seq) == 1 and
                low_count / high_count < percent, merge low to high.
        Returns:
            n_corrected_umi: int
            n_corrected_read: int
        """
        n_corrected_umi = 0
        n_corrected_read = 0

        # sort by value(UMI count) first, then key(UMI sequence)
        umi_arr = sorted(
            umi_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
        while True:
            # break when only highest in umi_arr
            if len(umi_arr) == 1:
                break
            umi_low = umi_arr.pop()
            low_seq = umi_low[0]
            low_count = umi_low[1]

            for umi_kv in umi_arr:
                high_seq = umi_kv[0]
                high_count = umi_kv[1]
                if float(low_count / high_count) > percent:
                    break
                if utils.hamming_distance(low_seq, high_seq) == 1:
                    n_low = umi_dict[low_seq]
                    n_corrected_umi += 1
                    n_corrected_read += n_low
                    # merge
                    umi_dict[high_seq] += n_low
                    del (umi_dict[low_seq])
                    break
        return n_corrected_umi, n_corrected_read

    @utils.add_log
    def bam2table(self):
        """
        bam to detail table
        must be used on name_sorted bam
        """
        samfile = pysam.AlignmentFile(self.bam, "rb")
        with open(self.count_detail_file, 'wt') as fh1:
            fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'count']) + '\n')

            def keyfunc(x):
                return x.query_name.split('_', 1)[0]
            for _, g in groupby(samfile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                for seg in g:
                    (barcode, umi) = seg.query_name.split('_')[:2]
                    if not seg.has_tag('XT'):
                        continue
                    gene_id = seg.get_tag('XT')
                    gene_umi_dict[gene_id][umi] += 1
                for gene_id in gene_umi_dict:
                    Count.correct_umi(gene_umi_dict[gene_id])

                # output
                for gene_id in gene_umi_dict:
                    for umi in gene_umi_dict[gene_id]:
                        fh1.write('%s\t%s\t%s\t%s\n' % (barcode, gene_id, umi,
                                                        gene_umi_dict[gene_id][umi]))
        samfile.close()

    @utils.add_log
    def cell_calling(self, df_sum):
        cell_calling_method = self.cell_calling_method

        if (self.force_cell_num is not None) and (self.force_cell_num != 'None'):
            cell_bc, UMI_threshold = self.force_cell(df_sum)
        elif cell_calling_method == 'auto':
            cell_bc, UMI_threshold = self.auto_cell(df_sum)
        elif cell_calling_method == 'EmptyDrops_CR':
            cell_bc, UMI_threshold = self.emptydrop_cr_cell(df_sum)
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
    def get_df_sum(df, col='UMI'):
        def num_gt2(x):
            return pd.Series.sum(x[x > 1])

        df_sum = df.groupby('Barcode').agg({
            'count': ['sum', num_gt2],
            'UMI': 'count',
            'geneID': 'nunique'
        })
        df_sum.columns = ['readcount', 'UMI2', 'UMI', 'geneID']
        df_sum = df_sum.sort_values(col, ascending=False)
        return df_sum

    def get_cell_stats(self, df_sum, cell_bc):
        df_sum.loc[:, 'mark'] = 'UB'
        df_sum.loc[df_sum.index.isin(cell_bc), 'mark'] = 'CB'
        df_sum.to_csv(self.marked_count_file, sep='\t')
        CB_describe = df_sum.loc[df_sum['mark'] == 'CB', :].describe()
        return CB_describe

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

    @utils.add_log
    def get_summary(self, CB_describe, CB_total_Genes,
                    CB_reads_count, reads_mapped_to_transcriptome):

        estimated_cells = int(CB_describe.loc['count', 'readcount'])
        self.add_metric(
            name='Estimated Number of Cells',
            value=estimated_cells,
            help_info=f'the number of barcodes considered as cell-associated. The cell calling method used for this sample is {self.cell_calling_method}.'
        )

        fraction_reads_in_cells = round(float(CB_reads_count) / reads_mapped_to_transcriptome * 100, 2)
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
            self.get_summary.logger.warning('Will not output `Mean Reads per Cell`')
        else:
            mean_reads_per_cell = int(valid_read_number / estimated_cells)
            self.add_metric(
                name='Mean Reads per Cell',
                value=mean_reads_per_cell,
                help_info='the number of valid reads divided by the estimated number of cells'
            )

        median_umi_per_cell = int(CB_describe.loc['50%', 'UMI'])
        self.add_metric(
            name='Median UMI per Cell',
            value=median_umi_per_cell,
            help_info='the median number of UMI counts per cell-associated barcode'
        )

        total_genes = int(CB_total_Genes)
        self.add_metric(
            name='Total Genes',
            value=total_genes,
            help_info='the number of genes with at least one UMI count in any cell'
        )

        median_genes_per_cell = int(CB_describe.loc['50%', 'geneID'])
        self.add_metric(
            name='Median Genes per Cell',
            value=median_genes_per_cell,
            help_info='the median number of genes detected per cell-associated barcode'
        )

        umi_saturation = round(self.downsample_dict['umi_saturation'][-1], 2)
        read_saturation = round(self.downsample_dict['read_saturation'][-1], 2)
        self.add_metric(
            name='Saturation',
            value=umi_saturation,
            display=f'{umi_saturation}%',
            help_info=(
                'the fraction of UMI originating from an already-observed UMI. '
                'There is a difference in how CeleScope and CellRanger calculate saturation. '
                'CeleScope shows umi_saturation in the report, while CellRanger shows read_saturation in the report. '
                'For details, see <a href="https://github.com/singleron-RD/CeleScope/blob/dev/docs/details.md#saturation">here</a>. '
                f'read_saturation: {read_saturation}%'
            )
        )

    @staticmethod
    def sub_sample(fraction, df_cell, cell_read_index):
        """
        umi_saturation = 1 - n_deduped_reads / n_umis
        read_saturation = 1 - n_deduped_reads / n_reads
        Currently the html report shows umi_saturation.

        n_deduped_reads = Number of unique (valid cell-barcode, valid UMI, gene) combinations among confidently mapped reads.
        n_umis = Total number of (confidently mapped, valid cell-barcode, valid UMI) UMIs.
        n_reads = Total number of (confidently mapped, valid cell-barcode, valid UMI) reads.
        Args:
            fration: subsmaple fration
            df_cell: in cell df with (Barcode geneID UMI count) 
            cell_read_index: df_cell repeat index
        """
        cell_read = df_cell['count'].sum()
        frac_n_read = int(cell_read * fraction)
        subsample_read_index = cell_read_index[:frac_n_read]
        index_dedup, counts = np.unique(subsample_read_index, return_counts=True)
        n_count_once = np.sum(counts == 1)
        # total = UMI
        umi_total = len(index_dedup)
        umi_saturation = round((1 - n_count_once / umi_total) * 100, 2)
        read_total = frac_n_read
        read_saturation = round((1 - n_count_once / read_total) * 100, 2)

        # gene median
        df_cell_subsample = df_cell.loc[index_dedup, ]
        geneNum_median = float(df_cell_subsample.groupby(
            'Barcode').agg({'geneID': 'nunique'}).median())

        return umi_saturation, read_saturation, geneNum_median

    @utils.add_log
    def downsample(self, df_cell):
        """saturation and median gene
        """
        cell_read_index = np.array(df_cell.index.repeat(df_cell['count']), dtype='int32')
        np.random.shuffle(cell_read_index)

        downsample_dict = {
            READ_FRACTION: [0],
            UMI_SATURATION: [0],
            READ_SATURATION: [0],
            MEDIAN_GENE_NUMBER: [0],
        }

        for fraction in np.arange(0.1, 1.1, 0.1):
            umi_saturation, read_saturation, geneNum_median = Count.sub_sample(
                fraction, df_cell, cell_read_index)
            fraction = round(fraction,1)
            umi_saturation = round(umi_saturation, 2)
            read_saturation = round(read_saturation, 2)
            downsample_dict[READ_FRACTION].append(fraction)
            downsample_dict[UMI_SATURATION].append(umi_saturation)
            downsample_dict[READ_SATURATION].append(read_saturation)
            downsample_dict[MEDIAN_GENE_NUMBER].append(geneNum_median)
        
            self.add_metric(
                name=f'Read Fraction {fraction} read_saturation',
                value=read_saturation,
                show=False,
            )

            self.add_metric(
                name=f'Read Fraction {fraction} umi_saturation',
                value=umi_saturation,
                show=False,
            )

        df_downsample = pd.DataFrame(downsample_dict, columns=[READ_FRACTION, MEDIAN_GENE_NUMBER, UMI_SATURATION, READ_SATURATION])
        df_downsample.to_csv(self.downsample_file, index=False, sep='\t')
        self.downsample_dict = downsample_dict


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
        parser.add_argument('--bam', help='Required. BAM file from featureCounts.', required=True)
        parser.add_argument(
            '--force_cell_num',
            help='Default `None`. Force the cell number to be this number. ',
        )


class Count_test(unittest.TestCase):
    def test_correct_umi(self):
        dic = {
            "apple1": 2,
            "apple2": 30,
            "bears1": 5,
            "bears2": 10,
            "bears3": 100,
            "ccccc1": 20,
            "ccccc2": 199,
        }
        n_corrected_umi, n_corrected_read = Count.correct_umi(dic)
        dic_after_correct = {
            'ccccc1': 20,
            'apple2': 32,
            'bears3': 115,
            'ccccc2': 199,
        }
        self.assertEqual(dic, dic_after_correct)
        self.assertEqual(n_corrected_umi, 3)
        self.assertEqual(n_corrected_read, 2 + 5 + 10)


if __name__ == "__main__":
    unittest.main()
