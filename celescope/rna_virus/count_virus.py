import logging
import os
from collections import defaultdict

import pandas as pd
import pysam

from celescope.tools.utils import add_log
from celescope.tools.step import s_common


def genDict(dim=3):
    if dim == 1:
        return defaultdict(int)
    else:
        return defaultdict(lambda: genDict(dim - 1))


def sum_virus(validated_barcodes, virus_bam,
              out_read_count_file, out_umi_count_file):
    # process bam
    samfile = pysam.AlignmentFile(virus_bam, "rb")
    count_dic = genDict(dim=3)
    for read in samfile:
        tag = read.reference_name
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if barcode in validated_barcodes:
            count_dic[barcode][tag][umi] += 1

    # write dic to pandas df
    rows = []
    for barcode in count_dic:
        for tag in count_dic[barcode]:
            for umi in count_dic[barcode][tag]:
                rows.append([barcode, tag, umi, count_dic[barcode][tag][umi]])
    if len(rows) == 0:
        logging.warning("No cell virus UMI found!")

    df_read = df_read = pd.DataFrame(
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
def count_virus(args):

    # 检查和创建输出目录
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % (args.outdir))

    # read barcodes
    df_barcodes = pd.read_csv(args.barcode_file, header=None)
    validated_barcodes = list(df_barcodes.iloc[:, 0])

    # count virus
    out_read_count_file = args.outdir + "/" + args.sample + "_virus_read_count.tsv"
    out_umi_count_file = args.outdir + "/" + args.sample + "_virus_UMI_count.tsv"
    sum_virus(
        validated_barcodes,
        args.virus_bam,
        out_read_count_file,
        out_umi_count_file)


def get_opts_count_virus(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument('--virus_bam', required=True)
        parser.add_argument('--barcode_file', required=True)
