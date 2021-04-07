import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import pysam
import os
import glob
from collections import defaultdict
from celescope.tools.utils import *
from celescope.tools.report import reporter
from itertools import groupby
from celescope.tools.count import correct_umi
import gzip
import math


def genDict(dim=3):
    if dim == 1:
        return defaultdict(int)
    else:
        return defaultdict(lambda: genDict(dim - 1))


@add_log
def sum_virus(validated_barcodes, virus_bam,
              out_read_count_file, out_umi_count_file, min_query_length, threshold=0.8):
    # process bam
    samfile = pysam.AlignmentFile(virus_bam, "rb")
    count_dic = genDict(dim=3)
    ref_num = samfile.nreferences
    for read in samfile:
        query_length = read.infer_query_length()
        tag = read.reference_name
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if (barcode in validated_barcodes) and (query_length >= min_query_length):
            count_dic[barcode][tag][umi] += 1

    ref_umi_dict = {barcode: count_dic[barcode] for barcode, count_dic[barcode] in count_dic.items()
                    if len(count_dic[barcode]) >= math.ceil(ref_num * threshold)}

    # write dic to pandas df
    rows = []
    for barcode in ref_umi_dict:
        for tag in ref_umi_dict[barcode]:
            for umi in ref_umi_dict[barcode][tag]:
                rows.append([barcode, tag, umi, ref_umi_dict[barcode][tag][umi]])
    if len(rows) == 0:
        sum_virus.logger.warning("No cell virus UMI found!")

    df_read = pd.DataFrame(
        rows,
        columns=[
            "barcode",
            "tag",
            "UMI",
            "read_count"])
    df_read.to_csv(out_read_count_file, sep="\t", index=False)

    df_umi = df_read.groupby(["barcode", "tag"]).agg({"UMI": "count"})
    df_umi.to_csv(out_umi_count_file, sep="\t")


@add_log
def count_capture_virus(args):

    # 检查和创建输出目录
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % (args.outdir))

    # read barcodes
    match_cell_barcodes, _match_cell_number = read_barcode_file(args.match_dir)

    # count virus
    out_read_count_file = args.outdir + "/" + args.sample + "_virus_read_count.tsv"
    out_umi_count_file = args.outdir + "/" + args.sample + "_virus_UMI_count.tsv"
    sum_virus(
        match_cell_barcodes,
        args.virus_bam,
        out_read_count_file,
        out_umi_count_file,
        args.min_query_length)


def get_opts_count_capture_virus(parser, sub_program):
    parser.add_argument(
        '--match_dir',
        help='matched rna_virus directory',
        required=True)
    parser.add_argument("--min_query_length", help='minimum query length', default=35)
    if sub_program:
        s_common(parser)
        parser.add_argument('--virus_bam', required=True)
