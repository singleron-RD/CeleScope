#!/bin/env python
# coding=utf8

from celescope.tools.utils import format_number
from celescope.tools.report import reporter
from collections import defaultdict
import pysam
import gzip
import os
import pandas as pd
import logging
import numpy as np
import sys
import argparse
import matplotlib as mpl
import re
import json
import glob
mpl.use('Agg')
from matplotlib import pyplot as plt

parentDir = os.path.dirname(__file__)
logger1 = logging.getLogger(__name__)


def genDict(dim=3):
    if dim == 1:
        return defaultdict(int)
    else:
        return defaultdict(lambda: genDict(dim - 1))


def read_pos(fusion_pos_file):
    df = pd.read_csv(fusion_pos_file, sep="\t")
    dic = dict([(tag, int(pos))
                for tag, pos in zip(df.iloc[:, 0], df.iloc[:, 1])])
    return dic


def is_fusion(pos, read_start, read_length, flanking_base):
    test_start = (pos - flanking_base) >= read_start
    test_end = (pos + flanking_base) <= (read_start + read_length)
    return (test_start and test_end)


def read_barcode_file(match_barcode_file):
    match_barcode = []
    with open(match_barcode_file, "rt") as f:
        while True:
            line = f.readline().strip()
            if not line:
                break
            match_barcode.append(line)
    return match_barcode


def count_fusion(args):

    logger1.info("count_fusion...!")

    outdir = args.outdir
    sample = args.sample
    bam = args.bam
    flanking_base = int(args.flanking_base)
    fusion_pos_file = args.fusion_pos
    match_dir = args.match_dir
    UMI_min = int(args.UMI_min)

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    fusion_pos = read_pos(fusion_pos_file)
    out_prefix = outdir + "/" + sample
    # barcode
    match_barcode_file1 = glob.glob(
        "{match_dir}/05.count/*_cellbarcode.tsv".format(match_dir=match_dir))
    match_barcode_file2 = glob.glob(
        "{match_dir}/05.count/matrix_10X/*_cellbarcode.tsv".format(match_dir=match_dir))
    match_barcode_file = (match_barcode_file1 + match_barcode_file2)[0]
    match_barcode = read_barcode_file(match_barcode_file)
    # tsne
    match_tsne_file = "{match_dir}/06.analysis/tsne_coord.tsv".format(
        match_dir=match_dir)
    df_tsne = pd.read_csv(match_tsne_file, sep="\t", index_col=0)
    # out
    out_read_count_file = out_prefix + "_fusion_read_count.tsv"
    out_umi_count_file = out_prefix + "_fusion_UMI_count.tsv"
    out_barcode_count_file = out_prefix + "_fusion_barcode_count.tsv"
    out_tsne_file = out_prefix + "_fusion_tsne.tsv"

    # process bam
    samfile = pysam.AlignmentFile(bam, "rb")
    header = samfile.header
    new_bam = pysam.AlignmentFile(
        out_prefix + "_fusion.bam", "wb", header=header)
    count_dic = genDict(dim=3)
    for read in samfile:
        tag = read.reference_name
        read_start = int(read.reference_start)
        read_length = len(read.query_sequence)
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if tag in fusion_pos.keys():
            if barcode in match_barcode:
                if is_fusion(pos=fusion_pos[tag], read_start=read_start,
                             read_length=read_length, flanking_base=flanking_base):
                    new_bam.write(read)
                    count_dic[barcode][tag][umi] += 1
    new_bam.close()

    # write dic to pandas df
    rows = []
    for barcode in count_dic:
        for tag in count_dic[barcode]:
            for umi in count_dic[barcode][tag]:
                rows.append([barcode, tag, umi, count_dic[barcode][tag][umi]])
    df_read = pd.DataFrame(rows)
    df_read.rename(
        columns={
            0: "barcode",
            1: "tag",
            2: "UMI",
            3: "read_count"},
        inplace=True)
    df_read.to_csv(out_read_count_file, sep="\t", index=False)

    df_umi = df_read.groupby(["barcode", "tag"]).agg({"UMI": "count"})
    df_umi = df_umi[df_umi["UMI"] >= UMI_min]
    df_umi.to_csv(out_umi_count_file, sep="\t")

    df_umi.reset_index(inplace=True)
    df_barcode = df_umi.groupby(["tag"]).agg({"barcode": "count"})
    n_match_barcode = len(match_barcode)
    # add zero count tag
    for tag in fusion_pos.keys():
        if not tag in df_barcode.barcode:
            new_row = pd.Series(data={'barcode': 0}, name=tag)
            df_barcode = df_barcode.append(new_row, ignore_index=False)
    df_barcode["percent"] = df_barcode["barcode"] / n_match_barcode
    df_barcode.to_csv(out_barcode_count_file, sep="\t")

    df_pivot = df_umi.pivot(index="barcode", columns="tag", values="UMI")
    df_pivot.fillna(0, inplace=True)
    df_tsne_fusion = pd.merge(
        df_tsne,
        df_pivot,
        right_index=True,
        left_index=True,
        how="left")
    df_tsne_fusion.fillna(0, inplace=True)
    df_tsne_fusion.to_csv(out_tsne_file, sep="\t")
    logger1.info("count done.")

    # plot
    logger1.info("plot fusion...!")
    app = parentDir + "/plot_fusion.R"
    cmd = f"Rscript {app} --tsne_fusion {out_tsne_file} --outdir {outdir}"
    os.system(cmd)
    logger1.info("plot done.")


def get_opts_count_fusion(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--bam", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    #parser.add_argument("--fusion_fasta",help="fusion fasta",required=True)
    parser.add_argument(
        "--fusion_pos",
        help="first base position of the second gene(0-start),tsv file",
        required=True)
    parser.add_argument(
        "--match_dir",
        help="match scRNA-Seq dir",
        required=True)
    parser.add_argument("--flanking_base", default=5)
    parser.add_argument("--UMI_min", default=1)
