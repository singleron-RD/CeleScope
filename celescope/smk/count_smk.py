import gzip
import os
import pandas as pd
import logging
import numpy as np
import sys
import argparse
import re
import glob
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from celescope.tools.utils import format_number
from celescope.tools.report import reporter

logger1 = logging.getLogger(__name__)


def tag_type(sub_df, UMI_min, percent_min):
    UMI_sum = sub_df["UMI_count"].sum()
    UMI_max = sub_df["UMI_count"].max()
    UMI_max_percent = UMI_max/UMI_sum
    if UMI_sum < UMI_min:
        return "Undetermined"
    elif UMI_max_percent < percent_min:
        return "Multiplet"
    else:
        return sub_df.loc[sub_df["UMI_count"] == UMI_max, "SMK_barcode"].values[0]


def write_and_plot(df, column_name, count_file, plot_file):
    df_count = df.groupby(["tag",column_name]).size().unstack()
    df_count.fillna(0, inplace=True)
    df_count.to_csv(count_file, sep="\t")
    df_percent = df_count/df_count.sum()
    df_plot = df_percent.stack().reset_index()
    df_plot.rename({0:"percent"}, axis=1, inplace=True)

    #plot
    colors = list(mpl.colors.cnames.keys())
    fig, ax = plt.subplots(figsize=(20, 10))  
    types = df_plot["tag"].drop_duplicates()
    margin_bottom = np.zeros(len(df_plot[column_name].drop_duplicates()))

    for num, type in enumerate(types):
        values = list(df_plot.loc[df_plot["tag"]==type, "percent"])
        df_plot[df_plot['tag'] == type].plot.bar(x=column_name, y='percent', ax=ax, stacked=True, 
            bottom=margin_bottom, label=type,color=colors[num*3+1])
        margin_bottom += values
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("SMK tag fraction")
    fig.savefig(plot_file)


def count_smk(args):

    logger1.info('count smk start...!')
    cell_UMI_file = args.cell_UMI_file
    tsne_file = glob.glob(f'{args.match_dir}/*analysis/tsne_coord.tsv')[0]
    UMI_min = int(args.UMI_min)
    percent_min = float(args.percent_min)
    combine_cluster = args.combine_cluster
    outdir = args.outdir
    sample = args.sample

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    tag_summary_file = f'{outdir}/{sample}_tag_summary.tsv'
    cluster_count_file = f'{outdir}/{sample}_cluster_count.tsv'
    cluster_plot = f'{outdir}/{sample}_cluster_plot.pdf'
    if combine_cluster:
        combine_cluster_count_file = f'{outdir}/{sample}_combine_cluster_count.tsv'
        combine_cluster_plot = f'{outdir}/{sample}_combine_cluster_plot.pdf'

    df_cell_UMI = pd.read_csv(cell_UMI_file, sep="\t", index_col=0)
    df_cell_UMI.fillna(0, inplace=True)

    df_tsne = pd.read_csv(tsne_file, sep="\t", index_col=0)

    df_tag_UMI = df_cell_UMI.stack().reset_index()
    df_tag_UMI.columns = ["barcode", "SMK_barcode", "UMI_count"]

    barcode_tag_type = df_tag_UMI.groupby(["barcode"]).apply(tag_type, UMI_min=UMI_min, 
        percent_min=percent_min)
    df_summary = pd.DataFrame(barcode_tag_type)
    df_summary.columns = ["tag"]
    df_tsne_summary = pd.merge(df_tsne, df_summary, how="left", left_index=True, right_index=True)
    df_tsne_summary.fillna({"tag": "Undetermined"}, inplace=True)
    if combine_cluster:
        df_combine_cluster = pd.read_csv(combine_cluster, sep="\t", header=None)
        df_combine_cluster.columns = ["cluster", "combine_cluster"]
        df_tsne_combine_cluster_summary = pd.merge(df_tsne_summary, df_combine_cluster,
            on=["cluster"], how="left", left_index=True).set_index(df_tsne_summary.index)
        df_tsne_combine_cluster_summary.to_csv(tag_summary_file, sep="\t")
    else:
        df_tsne_summary.to_csv(tag_summary_file, sep="\t")

    write_and_plot(df=df_tsne_summary, column_name="cluster", count_file=cluster_count_file,
        plot_file=cluster_plot)
    
    if combine_cluster:
        write_and_plot(df=df_tsne_combine_cluster_summary, column_name="combine_cluster",
        count_file=combine_cluster_count_file, plot_file=combine_cluster_plot)
    
    logger1.info('count smk done...!')


def get_opts_count_smk(parser, sub_program):
    parser.add_argument("--UMI_min", help="keep tags have UMI>=UMI_min", default = 50)
    parser.add_argument("--percent_min", help="keep tags have percent>=percent_min", default = 0.7)
    parser.add_argument("--combine_cluster", help="conbine cluster tsv file", default = None)
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path")
        parser.add_argument("--cell_UMI_file", help="cell SMK UMI file")