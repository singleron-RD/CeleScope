#!/bin/env python
# coding=utf8

import os
import sys
import json
import logging
import re
import numpy as np
import pandas as pd
import glob
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from celescope.tools.report import reporter
from celescope.tools.utils import glob_genomeDir, log, gene_convert
from celescope.rna.__init__ import __ASSAY__
from celescope.tools.Analysis import Analysis

toolsdir = os.path.dirname(__file__)


def report_prepare(outdir, tsne_df, marker_df):
    json_file = outdir + '/../.data.json'
    if not os.path.exists(json_file):
        data = {}
    else:
        fh = open(json_file)
        data = json.load(fh)
        fh.close()

    data["cluster_tsne"] = cluster_tsne_list(tsne_df)
    data["gene_tsne"] = gene_tsne_list(tsne_df)
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


def gene_tsne_list(tsne_df):
    """
    return data dic
    """
    tSNE_1 = list(tsne_df.tSNE_1)
    tSNE_2 = list(tsne_df.tSNE_2)
    Gene_Counts = list(tsne_df.Gene_Counts)
    res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "Gene_Counts": Gene_Counts}
    return res


def marker_table(marker_df):
    """
    return html code
    """
    marker_df = marker_df.loc[:, ["cluster", "gene",
                                  "avg_logFC", "pct.1", "pct.2", "p_val_adj"]]
    marker_df["cluster"] = marker_df["cluster"].apply(lambda x: f"cluster {x}")
    marker_gene_table = marker_df.to_html(
        escape=False,
        index=False,
        table_id="marker_gene_table",
        justify="center")
    return marker_gene_table


@log
def generate_matrix(gtf_file, matrix_file):

    id_name = gene_convert(gtf_file)
    matrix = pd.read_csv(matrix_file, sep="\t")

    gene_name_col = matrix.geneID.apply(lambda x: id_name[x])
    matrix.geneID = gene_name_col
    matrix = matrix.drop_duplicates(subset=["geneID"], keep="first")
    matrix = matrix.dropna()
    matrix = matrix.rename({"geneID": ""}, axis='columns')
    return matrix


@log
def seurat(sample, outdir, matrix_file, save_rds):
    app = toolsdir + "/run_analysis.R"
    cmd = (
        f'Rscript {app} --sample {sample} --outdir {outdir} --matrix_file {matrix_file} '
        f'--save_rds {save_rds}'
    )
    seurat.logger.info(cmd)
    os.system(cmd)


@log
def auto_assign(sample, outdir, type_marker_tsv):
    rds = f'{outdir}/{sample}.rds'
    app = toolsdir + "/auto_assign.R"
    cmd = (
        f'Rscript {app} '
        f'--rds {rds} '
        f'--type_marker_tsv {type_marker_tsv} '
        f'--outdir {outdir} '
        f'--sample {sample} '
    )
    auto_assign.logger.info(cmd)
    os.system(cmd)


@log
def analysis(args):

    # check dir
    outdir = args.outdir
    sample = args.sample
    assay = args.assay 

    matrix_file = args.matrix_file
    save_rds = args.save_rds
    type_marker_tsv = args.type_marker_tsv
    auto_assign_bool = False
    if type_marker_tsv and type_marker_tsv != 'None':
        auto_assign_bool = True
    if auto_assign_bool:
        save_rds = True

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # run_R
    seurat(sample, outdir, matrix_file, save_rds)

    # auto_assign
    if auto_assign_bool:
        auto_assign(sample, outdir, type_marker_tsv)

    # report
    ana = Analysis(     
        sample,
        outdir,
        assay,
        match_dir=outdir+'/../',
        step='analysis',     
    )
    ana.run()


def get_opts_analysis(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--matrix_file', help='matrix file', required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument('--save_rds', action='store_true', help='write rds to disk')
    parser.add_argument('--type_marker_tsv', help='cell type marker tsv')

