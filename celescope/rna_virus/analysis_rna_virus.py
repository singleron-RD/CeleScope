#!/bin/env python
# coding=utf8

import glob
import json
import os

import pandas as pd

import celescope.tools
from celescope.tools.report import reporter
from celescope.tools.utils import add_log, s_common

toolsdir = os.path.dirname(celescope.tools.__file__)


def report_prepare(outdir, tsne_df, marker_df, virus_df):
    json_file = outdir + '/../.data.json'
    if not os.path.exists(json_file):
        data = {}
    else:
        fh = open(json_file)
        data = json.load(fh)
        fh.close()

    data["cluster_tsne"] = cluster_tsne_list(tsne_df)
    data["virus_tsne"] = virus_tsne_list(tsne_df, virus_df)
    data["marker_gene_table"] = marker_table(marker_df)

    with open(json_file, 'w') as fh:
        json.dump(data, fh)


def cluster_tsne_list(tsne_df):
    """
    tSNE_1	tSNE_2	cluster Gene_Counts
    return data list
    """
    sum_df = tsne_df.groupby(["cluster"]).agg("count").iloc[:, 0]
    percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
    res = []
    for cluster in sorted(tsne_df.cluster.unique()):
        sub_df = tsne_df[tsne_df.cluster == cluster]
        name = "cluster {cluster}({percent}%)".format(
            cluster=cluster, percent=percent_df[cluster])
        tSNE_1 = list(sub_df.tSNE_1)
        tSNE_2 = list(sub_df.tSNE_2)
        res.append({"name": name, "tSNE_1": tSNE_1, "tSNE_2": tSNE_2})
    return res


def virus_tsne_list(tsne_df, virus_df):
    """
    return data dic
    """
    tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
    df = pd.merge(tsne_df, virus_df, on="barcode", how="left")
    df["UMI"] = df["UMI"].fillna(0)
    tSNE_1 = list(df.tSNE_1)
    tSNE_2 = list(df.tSNE_2)
    virus_UMI = list(df.UMI)
    res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "virus_UMI": virus_UMI}
    return res


def marker_table(marker_df):
    """
    return html code
    """
    marker_df = marker_df.loc[:, ["cluster", "gene",
                                  "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]]
    marker_gene_table = marker_df.to_html(
        escape=False,
        index=False,
        table_id="marker_gene_table",
        justify="center")
    return marker_gene_table


@add_log
def seurat(sample, outdir, matrix_file):
    app = toolsdir + "/run_analysis.R"
    cmd = f"Rscript {app} --sample {sample} --outdir {outdir} --matrix_file {matrix_file}"
    os.system(cmd)


@add_log
def analysis_rna_virus(args):

    # check dir
    outdir = args.outdir
    sample = args.sample
    matrix_file = args.matrix_file
    virus_file = args.virus_file

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # run_R
    seurat(sample, outdir, matrix_file)

    # report
    tsne_df_file = glob.glob("{outdir}/*tsne_coord.tsv".format(outdir=outdir))[0]
    marker_df_file = glob.glob("{outdir}/*markers.tsv".format(outdir=outdir))[0]
    tsne_df = pd.read_csv(tsne_df_file, sep="\t")
    marker_df = pd.read_csv(marker_df_file, sep="\t")
    virus_df = pd.read_csv(virus_file, sep="\t")

    report_prepare(outdir, tsne_df, marker_df, virus_df)

    t = reporter(
        name='analysis_rna_virus',
        assay=args.assay,
        sample=args.sample,
        outdir=args.outdir + '/..')
    t.get_report()


def get_opts_analysis_rna_virus(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument('--matrix_file', help='matrix file', required=True)
        parser.add_argument(
            '--virus_file',
            help='virus UMI count file',
            required=True)
        