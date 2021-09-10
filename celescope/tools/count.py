"""
count step
"""

import os
import random
import subprocess
import sys
from collections import defaultdict
from itertools import groupby

import numpy as np
import pandas as pd
import pysam
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

import celescope.tools.utils as utils
from celescope.tools.__init__ import (BARCODE_FILE_NAME, FEATURE_FILE_NAME,
                                      MATRIX_FILE_NAME)
from celescope.tools.cellranger3 import get_plot_elements
from celescope.tools.cellranger3.cell_calling_3 import cell_calling_3
from celescope.tools.step import Step, s_common
from celescope.rna.mkref import parse_genomeDir_rna

TOOLS_DIR = os.path.dirname(__file__)
random.seed(0)
np.random.seed(0)


class Count(Step):
    """
    Features
    - Cell-calling: Distinguish cell barcodes from background barcodes. 

    - Generate expression matrix.

    Output
    - `{sample}_all_matrix` The expression matrix of all detected barcodes. 
        Can be read in by calling the `Seurat::Read10X` function.

    - `{sample}_matrix_10X` The expression matrix of the barcode that is identified to be the cell. 
    Can be read in by calling the `Seurat::Read10X` function.

    - `{sample}_matrix.tsv.gz` The expression matrix of the barcode that is identified to be the cell, separated by tabs. 
    CeleScope >=1.2.0 does not output this file.

    - `{sample}_count_detail.txt.gz` 4 columns: 
        - barcode  
        - gene ID  
        - UMI count  
        - read_count  

    - `{sample}_counts.txt` 6 columns:
        - Barcode: barcode sequence
        - readcount: read count of each barcode
        - UMI2: UMI count (with reads per UMI >= 2) for each barcode
        - UMI: UMI count for each barcode
        - geneID: gene count for each barcode
        - mark: cell barcode or backgound barcode.

            `CB` cell  
            `UB` background  

    - `{sample}_downsample.txt` 3 columns：
        - percent: percentage of sampled reads
        - median_geneNum: median gene number per cell
        - saturation: sequencing saturation

    - `barcode_filter_magnitude.pdf` Barcode-UMI plot.


    """

    def __init__(self, args, step):
        Step.__init__(self, args, step)
        self.force_cell_num = args.force_cell_num
        self.cell_calling_method = args.cell_calling_method
        self.expected_cell_num = int(args.expected_cell_num)
        self.bam = args.bam

        # set
        self.gtf_file = parse_genomeDir_rna(args.genomeDir)['gtf']
        self.id_name = utils.get_id_name_dict(self.gtf_file)

        # output files
        self.count_detail_file = f'{self.outdir}/{self.sample}_count_detail.txt'
        self.marked_count_file = f'{self.outdir}/{self.sample}_counts.txt'
        self.raw_matrix_10X_dir = f'{self.outdir}/{self.sample}_all_matrix'
        self.cell_matrix_10X_dir = f'{self.outdir}/{self.sample}_matrix_10X'
        self.downsample_file = f'{self.outdir}/{self.sample}_downsample.txt'

    def run(self):
        self.bam2table()
        df = pd.read_table(self.count_detail_file, header=0)

        # df_sum
        df_sum = Count.get_df_sum(df)

        # export all matrix
        self.write_matrix_10X(df, self.raw_matrix_10X_dir)

        # call cells
        cell_bc, _threshold = self.cell_calling(df_sum)

        # get cell stats
        CB_describe = self.get_cell_stats(df_sum, cell_bc)

        # export cell matrix
        df_cell = df.loc[df['Barcode'].isin(cell_bc), :]
        self.write_matrix_10X(df_cell, self.cell_matrix_10X_dir)
        (CB_total_Genes, CB_reads_count, reads_mapped_to_transcriptome) = self.cell_summary(
            df, cell_bc)

        # downsampling
        cell_bc = set(cell_bc)
        saturation, res_dict = self.downsample(df_cell)

        # summary
        self.get_summary(saturation, CB_describe, CB_total_Genes,
                         CB_reads_count, reads_mapped_to_transcriptome)

        self.report_prepare()

        self.add_content_item('metric', downsample_summary=res_dict)
        self.clean_up()

    def report_prepare(self):

        df0 = pd.read_table(self.downsample_file, header=0)
        self.add_data_item(percentile=df0['percent'].tolist())
        self.add_data_item(MedianGeneNum=df0['median_geneNum'].tolist())
        self.add_data_item(Saturation=df0['saturation'].tolist())
        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(self.marked_count_file))
        self.add_data_item(umi_summary=True)

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
        elif cell_calling_method == 'cellranger3':
            cell_bc, UMI_threshold = self.cellranger3_cell(df_sum)
        elif cell_calling_method == 'inflection':
            _cell_bc, UMI_threshold = self.auto_cell(df_sum)
            cell_bc, UMI_threshold = self.inflection_cell(df_sum, UMI_threshold)
        return cell_bc, UMI_threshold

    @utils.add_log
    def force_cell(self, df_sum):
        force_cell_num = int(self.force_cell_num)
        cell_range = int(force_cell_num * 0.1)
        cell_low = force_cell_num - cell_range
        cell_high = force_cell_num + cell_range

        df_barcode_count = df_sum.groupby(
            ['UMI']).size().reset_index(
            name='barcode_counts')
        sorted_df = df_barcode_count.sort_values("UMI", ascending=False)
        sorted_df["barcode_cumsum"] = sorted_df["barcode_counts"].cumsum()
        for i in range(sorted_df.shape[0]):
            if sorted_df.iloc[i, :]["barcode_cumsum"] >= cell_low:
                index_low = i - 1
                break
        for i in range(sorted_df.shape[0]):
            if sorted_df.iloc[i, :]["barcode_cumsum"] >= cell_high:
                index_high = i
                break
        df_sub = sorted_df.iloc[index_low:index_high + 1, :]
        threshold = df_sub.iloc[np.argmax(
            np.diff(df_sub["barcode_cumsum"])), :]["UMI"]
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
    def cellranger3_cell(self, df_sum):
        cell_bc, initial_cell_num = cell_calling_3(self.raw_matrix_10X_dir, self.expected_cell_num)
        threshold = Count.find_threshold(df_sum, initial_cell_num)
        return cell_bc, threshold

    @utils.add_log
    def inflection_cell(self, df_sum, threshold):
        app = f'{TOOLS_DIR}/rescue.R'
        cmd = (
            f'Rscript {app} '
            f'--matrix_dir {self.raw_matrix_10X_dir} '
            f'--outdir {self.outdir} '
            f'--sample {self.sample} '
            f'--threshold {threshold}'
        )
        Count.inflection_cell.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
        out_file = f'{self.outdir}/{self.sample}_rescue.tsv'
        df = pd.read_csv(out_file, sep='\t')
        inflection = int(df.loc[:, 'inflection'])
        threshold = inflection
        cell_bc = Count.get_cell_bc(df_sum, threshold)

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

    '''
    @utils.add_log
    def plot_barcode_UMI(df_sum, threshold, expected_cell_num, cell_num, outdir, sample, cell_calling_method, col='UMI'):
        out_plot = f'{outdir}/{sample}_barcode_UMI_plot.pdf'
        import matplotlib
        import matplotlib.pyplot as plt
        fig = plt.figure()
        plt.plot(df_sum['UMI'])
        plt.hlines(threshold, 0, cell_num, linestyle='dashed')
        plt.vlines(cell_num, 0, threshold, linestyle='dashed')
        plt.title('cell_calling_method: %s, expected_cell_num: %s\n %s threshold: %s, cell num: %s' %
                (cell_calling_method, expected_cell_num, col, threshold, cell_num))
        plt.loglog()
        plt.savefig(out_plot)
    '''

    def get_cell_stats(self, df_sum, cell_bc):
        df_sum.loc[:, 'mark'] = 'UB'
        df_sum.loc[df_sum.index.isin(cell_bc), 'mark'] = 'CB'
        df_sum.to_csv(self.marked_count_file, sep='\t')
        CB_describe = df_sum.loc[df_sum['mark'] == 'CB', :].describe()
        return CB_describe

    @utils.add_log
    def write_matrix_10X(self, df, matrix_dir):
        if not os.path.exists(matrix_dir):
            os.mkdir(matrix_dir)

        df_UMI = df.groupby(['geneID', 'Barcode']).agg({'UMI': 'count'})
        mtx = coo_matrix((df_UMI.UMI, (df_UMI.index.codes[0], df_UMI.index.codes[1])))
        gene_id = df_UMI.index.levels[0].to_series()
        # add gene symbol
        gene_name = gene_id.apply(lambda x: self.id_name[x])
        genes = pd.concat([gene_id, gene_name], axis=1)
        genes.columns = ['gene_id', 'gene_name']

        barcodes = df_UMI.index.levels[1].to_series()
        genes.to_csv(f'{matrix_dir}/{FEATURE_FILE_NAME}', index=False, sep='\t', header=False)
        barcodes.to_csv(f'{matrix_dir}/{BARCODE_FILE_NAME}', index=False, sep='\t', header=False)
        mmwrite(f'{matrix_dir}/{MATRIX_FILE_NAME}', mtx)

    @utils.add_log
    def cell_summary(self, df, cell_bc):

        df.loc[:, 'mark'] = 'UB'
        df.loc[df['Barcode'].isin(cell_bc), 'mark'] = 'CB'
        CB_total_Genes = df.loc[df['mark'] == 'CB', 'geneID'].nunique()
        CB_reads_count = df.loc[df['mark'] == 'CB', 'count'].sum()
        reads_mapped_to_transcriptome = df['count'].sum()
        return(CB_total_Genes, CB_reads_count, reads_mapped_to_transcriptome)

    def get_summary(self, saturation, CB_describe, CB_total_Genes,
                    CB_reads_count, reads_mapped_to_transcriptome):

        # total read
        str_number = self.content_dict['data']['barcode_summary'][1][1].split("(")[0]
        valid_read_number = int(str_number.replace(",", ""))

        summary = pd.Series([0, 0, 0, 0, 0, 0, 0],
                            index=[
                                'Estimated Number of Cells',
                                'Fraction Reads in Cells',
                                'Mean Reads per Cell',
                                'Median UMI per Cell',
                                'Total Genes',
                                'Median Genes per Cell',
                                'Saturation',
        ])

        # 细胞数
        summary['Estimated Number of Cells'] = int(
            CB_describe.loc['count', 'readcount'])
        summary['Fraction Reads in Cells'] = '%.2f%%' % (float(
            CB_reads_count) / reads_mapped_to_transcriptome * 100)
        summary['Mean Reads per Cell'] = int(
            valid_read_number /
            summary['Estimated Number of Cells'])
        summary['Median UMI per Cell'] = int(CB_describe.loc['50%', 'UMI'])
        summary['Total Genes'] = int(CB_total_Genes)
        summary['Median Genes per Cell'] = int(CB_describe.loc['50%', 'geneID'])
        summary['Saturation'] = '%.2f%%' % (saturation)
        # 测序饱和度，认定为细胞中的reads中UMI>2的reads比例
        need_format = [
            'Estimated Number of Cells',
            'Mean Reads per Cell',
            'Median UMI per Cell',
            'Total Genes',
            'Median Genes per Cell']
        for item in need_format:
            summary[item] = utils.format_number(summary[item])
        summary.to_csv(self.stat_file, header=False, sep=':')

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
        return fraction=1 saturation
        """
        cell_read_index = np.array(df_cell.index.repeat(df_cell['count']), dtype='int32')
        np.random.shuffle(cell_read_index)

        format_str = "%.2f\t%.2f\t%.2f\n"
        res_dict = {
            "fraction": [],
            "umi_saturation": [],
            "read_saturation": [],
            "median_gene": []
        }
        with open(self.downsample_file, 'w') as fh:
            fh.write('percent\tmedian_geneNum\tsaturation\n')
            fh.write(format_str % (0, 0, 0))
            for fraction in np.arange(0.1, 1.1, 0.1):
                umi_saturation, read_saturation, geneNum_median = Count.sub_sample(
                    fraction, df_cell, cell_read_index)
                fh.write(format_str % (fraction, geneNum_median, umi_saturation))
                def format_float(x): return round(x / 100, 4)
                res_dict["fraction"].append(round(fraction, 1))
                res_dict["umi_saturation"].append(format_float(umi_saturation))
                res_dict["read_saturation"].append(format_float(read_saturation))
                res_dict["median_gene"].append(geneNum_median)

        return umi_saturation, res_dict


@utils.add_log
def count(args):
    step_name = "count"
    runner = Count(args, step_name)
    runner.run()


def get_opts_count(parser, sub_program):
    parser.add_argument('--genomeDir', help='Required. Genome directory.')
    parser.add_argument('--expected_cell_num', help='Default `3000`. Expected cell number.', default=3000)
    parser.add_argument(
        '--cell_calling_method',
        help='Default `auto`. Cell calling methods. Choose from `auto`, `cellranger3` and `inflection`.',
        choices=['auto', 'cellranger3', 'inflection', ],
        default='auto',
    )
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--bam', help='Required. BAM file from featureCounts.', required=True)
        parser.add_argument(
            '--force_cell_num',
            help='Default `None`. Force the cell number within (value * 0.9, value * 1.1). ',
        )
