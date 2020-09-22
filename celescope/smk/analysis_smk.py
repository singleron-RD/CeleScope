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
from celescope.tools.utils import glob_genomeDir

logger1 = logging.getLogger(__name__)


def report_prepare(outdir, tsne_df, marker_df, feature_df):
    json_file = outdir + '/../.data.json'
    if not os.path.exists(json_file):
        data = {}
    else:
        fh = open(json_file)
        data = json.load(fh)
        fh.close()

    data["cluster_tsne"] = cluster_tsne_list(tsne_df, "cluster")
    data["cluster_tsne_1"] = cluster_tsne_list(
        feature_df, "tag", show_tag=False)
    data["marker_gene_table"] = marker_table(marker_df)

    with open(json_file, 'w') as fh:
        json.dump(data, fh)


def cluster_tsne_list(tsne_df, colname, show_tag=True):
    """
    tSNE_1	tSNE_2	cluster Gene_Counts
    return data list
    """
    sum_df = tsne_df.groupby([colname]).agg("count").iloc[:, 0]
    percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
    res = []
    for cluster in sorted(tsne_df[colname].unique()):
        sub_df = tsne_df[tsne_df[colname] == cluster]
        if show_tag:
            name = f"{colname} {cluster}({percent_df[cluster]}%)"
        else:
            name = f"{cluster}({percent_df[cluster]}%)"
        tSNE_1 = list(sub_df.tSNE_1)
        tSNE_2 = list(sub_df.tSNE_2)
        res.append({"name": name, "tSNE_1": tSNE_1, "tSNE_2": tSNE_2})
    return res


def feature_tsne_list(feature_df):
    """
    return data dic
    """
    colname = "tag"
    tSNE_1 = list(feature_df.tSNE_1)
    tSNE_2 = list(feature_df.tSNE_2)
    feature_val = list(feature_df[colname])
    res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, colname: feature_val}
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


def analysis_smk(args):
    logger1.info('smk analysis ...!')

    # check dir
    outdir = args.outdir
    sample = args.sample
    tsne_tag_file = args.tsne_tag_file
    match_dir = args.match_dir

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # report
    tsne_df_file = glob.glob(f'{match_dir}/*analysis*/*tsne_coord.tsv')[0]
    marker_df_file = glob.glob(f'{match_dir}/*analysis*/*markers.tsv')[0]
    tsne_df = pd.read_csv(tsne_df_file, sep="\t")
    marker_df = pd.read_csv(marker_df_file, sep="\t")
    tsne_tag_df = pd.read_csv(tsne_tag_file, sep="\t", index_col=0)

    report_prepare(outdir, tsne_df, marker_df, tsne_tag_df)

    t = reporter(
        name='analysis_smk',
        assay=args.assay,
        sample=sample,
        outdir=args.outdir +
        '/..')
    t.get_report()


def get_opts_analysis_smk(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--match_dir', help='match_dir', required=True)
        parser.add_argument(
            '--tsne_tag_file',
            help='smk tsne tag file',
            required=True)
        parser.add_argument('--assay', help='assay', required=True)
