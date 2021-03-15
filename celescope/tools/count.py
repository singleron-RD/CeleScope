#!/bin/env python
# coding=utf8

import os
import sys
import json
import functools
import gzip
from collections import defaultdict
from itertools import groupby
import numpy as np
import pandas as pd
import subprocess
from scipy.io import mmwrite
from scipy.sparse import csr_matrix, coo_matrix
import pysam
from celescope.tools.utils import format_number, log, gene_convert, glob_genomeDir
from celescope.tools.report import reporter
from celescope.tools.cellranger3.cell_calling_3 import cell_calling_3
from celescope.tools.__init__ import MATRIX_FILE_NAME, FEATURE_FILE_NAME, BARCODE_FILE_NAME


toolsdir = os.path.dirname(__file__)


def report_prepare(count_file, downsample_file, outdir):

    json_file = outdir + '/.data.json'
    if not os.path.exists(json_file):
        data = {}
    else:
        fh = open(json_file)
        data = json.load(fh)
        fh.close()

    df0 = pd.read_table(downsample_file, header=0)
    data['percentile'] = df0['percent'].tolist()
    data['MedianGeneNum'] = df0['median_geneNum'].tolist()
    data['Saturation'] = df0['saturation'].tolist()

    #data['count' + '_summary'] = df0.T.values.tolist()

    df = pd.read_table(count_file, header=0)
    df = df.sort_values('UMI', ascending=False)
    data['CB_num'] = df[df['mark'] == 'CB'].shape[0]
    data['Cells'] = list(df.loc[df['mark'] == 'CB', 'UMI'])

    data['UB_num'] = df[df['mark'] == 'UB'].shape[0]
    data['Background'] = list(df.loc[df['mark'] == 'UB', 'UMI'])

    data['umi_summary'] = True

    with open(json_file, 'w') as fh:
        json.dump(data, fh)



def hd(x, y):
    return len([i for i in range(len(x)) if x[i] != y[i]])


def correct_umi(fh1, barcode, gene_umi_dict, percent=0.1):
    res_dict = defaultdict()

    for geneID in gene_umi_dict:
        _dict = gene_umi_dict[geneID]
        umi_arr = sorted(
            _dict.keys(), key=lambda x: (_dict[x], x), reverse=True)
        while True:
            # break when only one barcode or umi_low/umi_high great than 0.1
            if len(umi_arr) == 1:
                break
            umi_low = umi_arr.pop()

            for u in umi_arr:
                if float(_dict[umi_low]) / _dict[u] > percent:
                    break
                if hd(umi_low, u) == 1:
                    _dict[u] += _dict[umi_low]
                    del (_dict[umi_low])
                    break
        res_dict[geneID] = _dict
    return res_dict


@log
def bam2table(bam, detail_file):
    # 提取bam中相同barcode的reads，统计比对到基因的reads信息
    #
    samfile = pysam.AlignmentFile(bam, "rb")
    with gzip.open(detail_file, 'wt') as fh1:
        fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'count']) + '\n')

        # pysam.libcalignedsegment.AlignedSegment
        # AAACAGGCCAGCGTTAACACGACC_CCTAACGT_A00129:340:HHH72DSXX:2:1353:23276:30843
        # 获取read的barcode
        def keyfunc(x): return x.query_name.split('_', 1)[0]

        for _, g in groupby(samfile, keyfunc):
            gene_umi_dict = defaultdict(lambda: defaultdict(int))
            for seg in g:
                (barcode, umi) = seg.query_name.split('_')[:2]
                if not seg.has_tag('XT'):
                    continue
                geneID = seg.get_tag('XT')
                gene_umi_dict[geneID][umi] += 1
            res_dict = correct_umi(fh1, barcode, gene_umi_dict)

            # output
            for geneID in res_dict:
                for umi in res_dict[geneID]:
                    fh1.write('%s\t%s\t%s\t%s\n' % (barcode, geneID, umi,
                                                    res_dict[geneID][umi]))



@log
def cell_calling(cell_calling_method, force_cell_num, expected_cell_num, all_matrix_10X_dir, df_sum, outdir, sample):
    if (force_cell_num is not None) and (force_cell_num != 'None'):
        cell_bc, UMI_threshold = force_cell(force_cell_num, df_sum)
    elif cell_calling_method == 'auto':
        cell_bc, UMI_threshold = auto_cell(df_sum, expected_cell_num)
    elif cell_calling_method == 'cellranger3':
        cell_bc, UMI_threshold = cellranger3_cell(all_matrix_10X_dir, expected_cell_num, df_sum)
    elif cell_calling_method == 'inflection':
        _cell_bc, UMI_threshold = auto_cell(df_sum, expected_cell_num)
        cell_bc, UMI_threshold = inflection_cell(outdir, sample, all_matrix_10X_dir, df_sum, UMI_threshold)
    cell_calling.logger.info(f'UMI_threshold: {UMI_threshold}')
    return cell_bc, UMI_threshold


@log
def force_cell(force_cell_num, df_sum):
    force_cell_num = int(force_cell_num)
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
    cell_bc = get_cell_bc(df_sum, threshold, col='UMI')

    return cell_bc, threshold


def find_threshold(df_sum, idx):
    return int(df_sum.iloc[idx - 1, df_sum.columns == 'UMI'])


@log
def auto_cell(df_sum, expected_cell_num):
    col = "UMI"
    idx = int(expected_cell_num * 0.01)
    barcode_number = df_sum.shape[0]
    idx = int(min(barcode_number, idx))
    if idx == 0:
        sys.exit("cell number equals zero!")
    # calculate read counts threshold
    threshold = int(find_threshold(df_sum, idx) * 0.1)
    threshold = max(1, threshold)
    cell_bc = get_cell_bc(df_sum, threshold)

    return cell_bc, threshold


@log
def cellranger3_cell(all_matrix_10X_dir, expected_cell_num, df_sum):
    cell_bc, initial_cell_num = cell_calling_3(all_matrix_10X_dir, expected_cell_num)
    threshold = find_threshold(df_sum, initial_cell_num)
    return cell_bc, threshold


@log
def inflection_cell(outdir, sample, all_matrix_10X_dir, df_sum, threshold):
    app = f'{toolsdir}/rescue.R'
    cmd = (
        f'Rscript {app} '
        f'--matrix_dir {all_matrix_10X_dir} '
        f'--outdir {outdir} '
        f'--sample {sample} '
        f'--threshold {threshold}'
    )
    inflection_cell.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)
    out_file = f'{outdir}/{sample}_rescue.tsv'
    df = pd.read_csv(out_file, sep='\t')
    inflection = int(df.loc[:,'inflection'])
    threshold = inflection
    cell_bc = get_cell_bc(df_sum, threshold)

    return cell_bc, threshold


@log
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

def get_cell_bc(df_sum, threshold, col='UMI'):
    return list(df_sum[df_sum[col] >= threshold].index)

@log
def plot_barcode_UMI(df_sum, threshold, expected_cell_num, cell_num, outdir, sample, cell_calling_method, col='UMI'):
    out_plot = f'{outdir}/{sample}_barcode_UMI_plot.pdf'
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    plt.plot(df_sum['UMI'])
    plt.hlines(threshold, 0, cell_num, linestyle='dashed')
    plt.vlines(cell_num, 0, threshold, linestyle='dashed')
    plt.title('cell_calling_method: %s, expected_cell_num: %s\n %s threshold: %s, cell num: %s' %
              (cell_calling_method, expected_cell_num, col, threshold, cell_num))
    plt.loglog()
    plt.savefig(out_plot)


def get_cell_stats(df_sum, cell_bc, marked_counts_file):
    df_sum.loc[:, 'mark'] = 'UB'
    df_sum.loc[df_sum.index.isin(cell_bc), 'mark'] = 'CB'
    df_sum.to_csv(marked_counts_file, sep='\t')
    CB_describe = df_sum.loc[df_sum['mark'] == 'CB', :].describe()

    return CB_describe


def write_matrix_10X(table, id_name, matrix_10X_dir):
    id = table.index.to_series()
    name = id.apply(lambda x: id_name[x])
    genes = pd.concat([id, name], axis=1)
    genes.columns = ['gene_id', 'gene_name']

    #write
    table.columns.to_series().to_csv(
        f'{matrix_10X_dir}/barcodes.tsv', index=False, sep='\t')
    genes.to_csv(
        f'{matrix_10X_dir}/genes.tsv', index=False, header=False, sep='\t')
    mmwrite(f'{matrix_10X_dir}/matrix', csr_matrix(table))
    return id, name


@log
def matrix_10X(df, outdir, sample, gtf_file, dir_name='matrix_10X', cell_bc=None):
    matrix_10X_dir = f"{outdir}/{sample}_{dir_name}/"
    if not os.path.exists(matrix_10X_dir):
        os.mkdir(matrix_10X_dir)
    id_name = gene_convert(gtf_file)

    if cell_bc is not None:
        df = df.loc[df['Barcode'].isin(cell_bc), :]
    
    df_UMI = df.groupby(['geneID','Barcode']).agg({'UMI':'count'})
    mtx= coo_matrix((df_UMI.UMI, (df_UMI.index.labels[0], df_UMI.index.labels[1])))
    id = df_UMI.index.levels[0].to_series()
    # add gene symbol
    name = id.apply(lambda x: id_name[x])
    genes = pd.concat([id, name], axis=1)
    genes.columns = ['gene_id', 'gene_name']

    barcodes = df_UMI.index.levels[1].to_series()
    genes.to_csv(f'{matrix_10X_dir}/{FEATURE_FILE_NAME}', index=False, sep='\t', header=False)
    barcodes.to_csv(f'{matrix_10X_dir}/{BARCODE_FILE_NAME}', index=False, sep='\t')
    mmwrite(f'{matrix_10X_dir}/{MATRIX_FILE_NAME}', mtx)
    return matrix_10X_dir


@log
def expression_matrix(df, cell_bc, outdir, sample, gtf_file):

    id_name = gene_convert(gtf_file)

    df.loc[:, 'mark'] = 'UB'
    df.loc[df['Barcode'].isin(cell_bc), 'mark'] = 'CB'
    CB_total_Genes = df.loc[df['mark'] == 'CB', 'geneID'].nunique()
    CB_reads_count = df.loc[df['mark'] == 'CB', 'count'].sum()
    reads_mapped_to_transcriptome = df['count'].sum()

    table = df.loc[df['mark'] == 'CB', :].pivot_table(
        index='geneID', columns='Barcode', values='UMI',
        aggfunc=len).fillna(0).astype(int)

    # convert id to name; write table matrix
    matrix_table_file = f"{outdir}/{sample}_matrix.tsv.gz"
    id = table.index.to_series()
    name = id.apply(lambda x: id_name[x])
    table.index = name
    table.index.name = ""
    table.to_csv(
        matrix_table_file,
        sep="\t",
        compression='gzip')
    return(CB_total_Genes, CB_reads_count, reads_mapped_to_transcriptome)


def get_summary(df, sample, Saturation, CB_describe, CB_total_Genes,
                CB_reads_count, reads_mapped_to_transcriptome,
                stat_file, outdir):

    # total read
    json_file = outdir + '.data.json'
    fh = open(json_file)
    data = json.load(fh)
    str_number = data['barcode_summary'][1][1].split("(")[0]
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
    summary['Saturation'] = '%.2f%%' % (Saturation)
    # 测序饱和度，认定为细胞中的reads中UMI>2的reads比例
    need_format = [
        'Estimated Number of Cells',
        'Mean Reads per Cell',
        'Median UMI per Cell',
        'Total Genes',
        'Median Genes per Cell']
    for item in need_format:
        summary[item] = format_number(summary[item])
    summary.to_csv(stat_file, header=False, sep=':')


def sample(p, df, bc):
    tmp = df.sample(frac=p)
    format_str = "%.2f\t%.2f\t%.2f\n"
    tmp.reset_index(drop=True, inplace=True)
    tmp = tmp.loc[tmp['Barcode'].isin(bc), :]
    geneNum_median = tmp.groupby('Barcode').agg({
        'geneID': 'nunique'
    })['geneID'].median()
    tmp = tmp.groupby(['Barcode', 'geneID', 'UMI']).agg({'UMI': 'count'})
    total = float(tmp['UMI'].count())
    saturation = (1 - tmp.loc[tmp['UMI'] == 1, :].shape[0] / total) * 100
    return format_str % (p, geneNum_median, saturation), saturation


@log
def downsample(detail_file, cell_bc, downsample_file):
    df = pd.read_table(detail_file, index_col=[0, 1, 2])
    df = df.index.repeat(df['count']).to_frame()
    format_str = "%.2f\t%.2f\t%.2f\n"
    with open(downsample_file, 'w') as fh:
        fh.write('percent\tmedian_geneNum\tsaturation\n')
        fh.write(format_str % (0, 0, 0))

        saturation = 0
        for p in np.arange(0.1, 1.1, 0.1):
            r, s = sample(p=p, df=df, bc=cell_bc)
            fh.write(r)
            saturation = s
    return saturation


@log
def count(args):
    # args
    outdir = args.outdir
    sample = args.sample
    assay = args.assay
    force_cell_num = args.force_cell_num
    cell_calling_method = args.cell_calling_method
    expected_cell_num = int(args.expected_cell_num)

    # check
    if args.genomeDir and args.genomeDir != "None":
        _refFlat, gtf_file = glob_genomeDir(args.genomeDir)
    else:
        gtf_file = args.gtf

    # 检查和创建输出目录
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # umi纠错，输出Barcode geneID  UMI     count为表头的表格
    count_detail_file = outdir + '/' + sample + '_count_detail.txt.gz'
    bam2table(args.bam, count_detail_file)

    df = pd.read_table(count_detail_file, header=0)

    # export all matrix
    dir_name = 'all_matrix'
    all_matrix_10X_dir = matrix_10X(df, outdir, sample, gtf_file, dir_name=dir_name)

    # df_sum
    df_sum = get_df_sum(df)

    # call cells
    cell_bc, threshold = cell_calling(cell_calling_method, force_cell_num, expected_cell_num, all_matrix_10X_dir, df_sum, outdir, sample)

    # plot
    cell_num = len(cell_bc)
    plot_barcode_UMI(df_sum, threshold, expected_cell_num, cell_num, outdir, sample, cell_calling_method, col='UMI')
    
    # get cell stats
    marked_counts_file = outdir + '/' + sample + '_counts.txt'
    CB_describe = get_cell_stats(df_sum, cell_bc, marked_counts_file)

    # export cell matrix
    matrix_10X(df, outdir, sample, gtf_file, dir_name='matrix_10X', cell_bc=cell_bc)
    (CB_total_Genes, CB_reads_count,
        reads_mapped_to_transcriptome) = expression_matrix(
            df, cell_bc, outdir, sample, gtf_file)

    # downsampling
    cell_bc = set(cell_bc)
    downsample_file = args.outdir + '/' + args.sample + '_downsample.txt'
    Saturation = downsample(
        count_detail_file,
        cell_bc,
        downsample_file)

    # summary
    stat_file = outdir + '/stat.txt'
    get_summary(df, sample, Saturation, CB_describe, CB_total_Genes,
                CB_reads_count, reads_mapped_to_transcriptome, stat_file,
                outdir + '/../')

    report_prepare(marked_counts_file, downsample_file, outdir + '/..')

    t = reporter(assay=assay,
                 name='count', sample=args.sample,
                 stat_file=outdir + '/stat.txt',
                 outdir=outdir + '/..')
    t.get_report()


def get_opts_count(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--bam', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument('--genomeDir', help='genome directory')
        parser.add_argument('--gtf', help='gtf file path')
    parser.add_argument('--force_cell_num', help='force cell number', default=None)
    parser.add_argument('--expected_cell_num', help='expected cell number', default=3000)
    parser.add_argument('--cell_calling_method', help='cell calling methods', 
        choices=['auto', 'cellranger3', 'inflection',], default='auto')