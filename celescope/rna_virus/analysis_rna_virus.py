#!/bin/env python
# coding=utf8

import os
import sys
import json
import logging
import re
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from celescope.tools.report import reporter
from celescope.tools.utils import glob_genomeDir

logger1 = logging.getLogger(__name__)
# invoke by celescope
toolsdir = os.path.dirname(__file__) + "../tools/"


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
                                  "avg_logFC", "pct.1", "pct.2", "p_val_adj"]]
    marker_gene_table = marker_df.to_html(
        escape=False,
        index=False,
        table_id="marker_gene_table",
        justify="center")
    return marker_gene_table


def gene_convert(gtf_file, matrix_file):

    gene_id_pattern = re.compile(r'gene_id "(\S+)";')
    gene_name_pattern = re.compile(r'gene_name "(\S+)"')
    id_name = {}
    with open(gtf_file) as f:
        for line in f.readlines():
            if line.startswith('#!'):
                continue
            tabs = line.split('\t')
            gtf_type, attributes = tabs[2], tabs[-1]
            if gtf_type == 'gene':
                gene_id = gene_id_pattern.findall(attributes)[-1]
                gene_name = gene_name_pattern.findall(attributes)[-1]
                id_name[gene_id] = gene_name

    matrix = pd.read_csv(matrix_file, sep="\t")

    def convert(gene_id):
        if gene_id in id_name:
            return id_name[gene_id]
        else:
            return np.nan
    gene_name_col = matrix.geneID.apply(convert)
    matrix.geneID = gene_name_col
    matrix = matrix.drop_duplicates(subset=["geneID"], keep="first")
    matrix = matrix.dropna()
    matrix = matrix.rename({"geneID": ""}, axis='columns')
    return matrix


def analysis_rna_virus(args):
    logger1.info('virus_analysis ...!')

    refFlat, gtf_file = glob_genomeDir(args.genomeDir, logger1)
    # check dir
    outdir = args.outdir
    sample = args.sample
    matrix_file = args.matrix_file
    virus_file = args.virus_file

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # run
    logger1.info("convert expression matrix.")
    new_matrix = gene_convert(gtf_file, matrix_file)
    new_matrix_file = "{outdir}/{sample}_matrix.tsv.gz".format(
        outdir=outdir, sample=sample)
    new_matrix.to_csv(
        new_matrix_file,
        sep="\t",
        index=False,
        compression='gzip')
    logger1.info("expression matrix written.")

    # run_R
    logger1.info("Seurat running.")
    cmd = "Rscript {app} --sample {sample} --outdir {outdir} --matrix_file {new_matrix_file}".format(
        app=toolsdir + "/run_analysis.R", sample=sample, outdir=outdir, new_matrix_file=new_matrix_file)
    os.system(cmd)
    logger1.info("Seurat done.")

    # report
    tsne_df_file = "{outdir}/tsne_coord.tsv".format(outdir=outdir)
    marker_df_file = "{outdir}/markers.tsv".format(outdir=outdir)
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
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument(
            '--matrix_file',
            help='matrix file xls',
            required=True)
        parser.add_argument(
            '--genomeDir',
            help='genome directory',
            required=True)
        parser.add_argument(
            '--virus_file',
            help='virus UMI count file',
            required=True)
        parser.add_argument('--assay', help='assay', required=True)
