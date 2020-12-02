from celescope.tools.report import reporter
from celescope.tools.utils import format_number, read_barcode_file, log, format_stat
import matplotlib.pyplot as plt
import os
import pandas as pd
import logging
import numpy as np
import argparse
import glob
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')


def get_UMI(row):
    return row.sum()


def get_UMI_min(df_cell_UMI, UMI_min):
    if UMI_min == "auto":
        UMI_min1 = np.percentile(df_cell_UMI.sum(axis=1), 5)
        UMI_min2 = np.median(df_cell_UMI.sum(axis=1)) / 10
        UMI_min = int(min(UMI_min1, UMI_min2))
        UMI_min = max(UMI_min, 1)
        return UMI_min
    else:
        return int(UMI_min)


def get_SNR(row, dim):
    row_sorted = sorted(row, reverse=True)
    noise = row_sorted[dim]
    signal = row_sorted[dim - 1]
    if noise == 0:
        return np.inf
    return float(signal) / noise


def get_SNR_min(df_cell_UMI, dim, SNR_min, UMI_min):
    UMIs = df_cell_UMI.apply(get_UMI, axis=1)
    df_valid_cell_UMI = df_cell_UMI[UMIs >= UMI_min]
    if SNR_min == "auto":
        SNRs = df_valid_cell_UMI.apply(get_SNR, dim=dim, axis=1)
        if np.median(SNRs) == np.inf:
            return 10
        return max(np.median(SNRs) / 10, 2)
    else:
        return float(SNR_min)


def tag_type(row, UMI_min, SNR_min, dim):
    SNR = get_SNR(row, dim)
    UMI = get_UMI(row)
    if UMI < UMI_min:
        return "Undetermined"
    if SNR < SNR_min:
        return "Multiplet"
    # get tag
    signal_tags = sorted(row.sort_values(ascending=False).index[0:dim])
    signal_tags_str = "_".join(signal_tags)
    return signal_tags_str


def write_and_plot(df, column_name, count_file, plot_file):
    df_count = df.groupby(["tag", column_name]).size().unstack()
    df_count.fillna(0, inplace=True)
    df_count.to_csv(count_file, sep="\t")
    df_percent = df_count / df_count.sum()
    df_plot = df_percent.stack().reset_index()
    df_plot.rename({0: "percent"}, axis=1, inplace=True)

    # plot
    colors = list(mpl.colors.cnames.keys())
    fig, ax = plt.subplots(figsize=(20, 10))
    types = df_plot["tag"].drop_duplicates()
    margin_bottom = np.zeros(len(df_plot[column_name].drop_duplicates()))

    for num, type in enumerate(types):
        values = list(df_plot.loc[df_plot["tag"] == type, "percent"])
        df_plot[df_plot['tag'] == type].plot.bar(
            x=column_name, y='percent', ax=ax, stacked=True,
            bottom=margin_bottom, label=type, color=colors[num * 3 + 1])
        margin_bottom += values
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("SMK tag fraction")
    fig.savefig(plot_file)


@log
def count_smk(args):

    read_file = args.read_file
    match_dir = args.match_dir
    tsne_file = glob.glob(f'{match_dir}/*analysis/*tsne_coord.tsv')[0]
    UMI_min = args.UMI_min
    SNR_min = args.SNR_min
    dim = int(args.dim)
    combine_cluster = args.combine_cluster
    outdir = args.outdir
    sample = args.sample
    assay = args.assay

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # stat_row
    stats = pd.Series()

    # process
    match_barcode, cell_total = read_barcode_file(match_dir)

    UMI_tag_file = f'{outdir}/{sample}_umi_tag.tsv'
    tsne_tag_file = f'{outdir}/{sample}_tsne_tag.tsv'
    cluster_count_file = f'{outdir}/{sample}_cluster_count.tsv'
    cluster_plot = f'{outdir}/{sample}_cluster_plot.pdf'
    if combine_cluster:
        combine_cluster_count_file = f'{outdir}/{sample}_combine_cluster_count.tsv'
        combine_cluster_plot = f'{outdir}/{sample}_combine_cluster_plot.pdf'

    df_read_count = pd.read_csv(read_file, sep="\t", index_col=0)
    mapped_read = df_read_count['read_count'].sum()

    # in cell
    df_read_count_in_cell = df_read_count[df_read_count.index.isin(match_barcode)]
    mapped_read_in_cell = int(df_read_count_in_cell['read_count'].sum())
    stats = stats.append(pd.Series(
        format_stat(mapped_read_in_cell, mapped_read),
        index=['Mapped Reads in Cells']
    ))

    # UMI
    df_UMI_in_cell = df_read_count_in_cell.reset_index().groupby([
        'barcode', 'SMK_barcode_name']).agg({'UMI': 'count'})
    df_UMI_in_cell = df_UMI_in_cell.reset_index()
    df_UMI_in_cell = df_UMI_in_cell.pivot(
        index='barcode', columns='SMK_barcode_name', values='UMI')
    df_cell = pd.DataFrame(index=match_barcode)
    df_UMI_cell = pd.merge(
        df_cell,
        df_UMI_in_cell,
        how="left",
        left_index=True,
        right_index=True)

    # fillna
    df_UMI_cell.fillna(0, inplace=True)
    df_UMI_cell = df_UMI_cell.astype(int)

    # UMI
    UMIs = df_UMI_cell.apply(sum, axis=1)
    median = round(np.median(UMIs), 2)
    mean = round(np.mean(UMIs), 2)
    stats = stats.append(pd.Series(
        str(median),
        index=['Median UMI per Cell']
    ))

    stats = stats.append(pd.Series(
        str(mean),
        index=['Mean UMI per Cell']
    ))

    UMI_min = get_UMI_min(df_UMI_cell, UMI_min)
    count_smk.logger.info(f'UMI_min: {UMI_min}')
    SNR_min = get_SNR_min(df_UMI_cell, dim, SNR_min, UMI_min)
    count_smk.logger.info(f'SNR_min: {SNR_min}')
    df_UMI_cell["tag"] = df_UMI_cell.apply(
        tag_type, UMI_min=UMI_min, SNR_min=SNR_min, dim=dim, axis=1)
    df_UMI_cell.to_csv(UMI_tag_file, sep="\t")

    df_tsne = pd.read_csv(tsne_file, sep="\t", index_col=0)
    df_tsne_tag = pd.merge(
        df_tsne,
        df_UMI_cell,
        how="left",
        left_index=True,
        right_index=True)

    if combine_cluster:
        df_combine_cluster = pd.read_csv(
            combine_cluster, sep="\t", header=None)
        df_combine_cluster.columns = ["cluster", "combine_cluster"]
        df_tsne_combine_cluster_tag = pd.merge(
            df_tsne_tag, df_combine_cluster,
            on=["cluster"], how="left", left_index=True).set_index(df_tsne_tag.index)
        df_tsne_combine_cluster_tag.to_csv(tsne_tag_file, sep="\t")
    else:
        df_tsne_tag.to_csv(tsne_tag_file, sep="\t")

    write_and_plot(
        df=df_tsne_tag, column_name="cluster", count_file=cluster_count_file,
        plot_file=cluster_plot
    )

    if combine_cluster:
        write_and_plot(
            df=df_tsne_combine_cluster_tag,
            column_name="combine_cluster",
            count_file=combine_cluster_count_file,
            plot_file=combine_cluster_plot
        )

    df_tag_count = df_UMI_cell["tag"].value_counts().reset_index()
    df_tag_count.columns = ["item", "count"]
    for index, row in df_tag_count.iterrows():
        stats = stats.append(pd.Series(
            format_stat(row['count'], cell_total),
            index=[row['item'] + ' Cells']
    ))
    stat_file = f'{outdir}/stat.txt'
    stats.to_csv(stat_file, sep=':', header=False)


    t = reporter(
        name='count_smk', assay=assay, sample=sample,
        stat_file=stat_file, outdir=outdir + '/..')
    t.get_report()


def get_opts_count_smk(parser, sub_program):
    parser.add_argument(
        "--UMI_min",
        help="cells have SMK_UMI>=UMI_min are considered as valid cell",
        default="auto")
    parser.add_argument("--dim", help="SMK tag dimension", default=1)
    parser.add_argument(
        "--SNR_min",
        help="minimum signal to noise ratio",
        default="auto")
    parser.add_argument(
        "--combine_cluster",
        help="conbine cluster tsv file",
        default=None)
    parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path", required=True)
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument("--read_file", help="SMK read count file")
