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
from celescope.tools.utils import log, format_number, genDict, parse_match_dir, read_barcode_file
import celescope.fusion

fusionDir = os.path.dirname(celescope.fusion.__file__)


def read_pos(fusion_pos_file):
    dic = {}
    df = pd.read_csv(fusion_pos_file, sep="\t")
    for tag, pos in zip(df.iloc[:, 0], df.iloc[:, 1]):
        dic[tag] = pos
    return dic


def is_fusion(pos, read_start, read_length, flanking_base):
    test_start = (pos - flanking_base) >= read_start
    test_end = (pos + flanking_base) <= (read_start + read_length)
    return (test_start and test_end)


@log
def count_fusion(args):

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
    match_barcode, _n_barcode = read_barcode_file(match_dir)
    # tsne
    match_tsne_file = parse_match_dir(match_dir)['tsne_coord']
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

    if not rows:
        count_fusion.logger.error('***** NO FUSION FOUND! *****')
    else:
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

        # plot
        count_fusion.logger.info("plot fusion...!")
        app = fusionDir + "/plot_fusion.R"
        cmd = f"Rscript {app} --tsne_fusion {out_tsne_file} --outdir {outdir}"
        os.system(cmd)
        count_fusion.logger.info("plot done.")


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
