#!/bin/env python
#coding=utf8

import os
import sys
import json
import functools
import logging
from collections import defaultdict
from itertools import groupby

import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import pysam

FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level = logging.INFO, format = FORMAT)

def report_prepare(outdir,tsne_df):
    json_file = outdir + '/.data.json'
    if not os.path.exists(json_file):
        data = {}
    else:
        fh = open(json_file)
        data = json.load(fh)
        fh.close()

    data["cluster_tsne"] = cluster_tsne_list(tsne_df)
    data["gene_tsne"] = gene_tsne_list(tsne_df)

    with open(json_file, 'w') as fh:
        json.dump(data, fh)

def cluster_tsne_list(tsne_df):
    """
    tSNE_1	tSNE_2	cluster Gene_Counts
    return data list
    """
    tsne_df.cluster = tsne_df.cluster + 1
    res = []
    for cluster in sorted(tsne_df.cluster.unique()):
        sub_df = tsne_df[tsne_df.cluster==cluster]
        name = "cluster" + str(cluster)
        tSNE_1 = list(sub_df.tSNE_1)
        tSNE_2 = list(sub_df.tSNE_2)
        res.append({"name":name,"tSNE_1":tSNE_1,"tSNE_2":tSNE_2})
    return res

def gene_tsne_list(tsne_df):
    """
    return data dic
    """
    tSNE_1 = list(tsne_df.tSNE_1)
    tSNE_2 = list(tsne_df.tSNE_2)
    Gene_Counts = list(tsne_df.Gene_Counts)
    res = {"tSNE_1":tSNE_1,"tSNE_2":tSNE_2,"Gene_Counts":Gene_Counts}
    return res

def marker_table(marker_df):
    pass



def analysis():
    pass


if __name__ == "__main__":
    tsne_df = pd.read_csv("/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/scope_tools_1.0/out/06.analysis/tsne_coord.tsv",sep="\t")
    report_prepare("./out",tsne_df)
    from report import reporter
    t = reporter(
        name='analysis',
        #stat_file=args.outdir + '/stat.txt',
        outdir= './out')
    t.get_report()
