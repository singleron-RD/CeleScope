import os

import pandas as pd
import pysam

from celescope.tools.utils import add_log, genDict, read_barcode_file, s_common


@add_log
def sum_virus(validated_barcodes, virus_bam,
              out_read_count_file, out_umi_count_file, min_query_length):
    # process bam
    samfile = pysam.AlignmentFile(virus_bam, "rb")
    count_dic = genDict(dim=3)
    for read in samfile:
        tag = read.reference_name
        query_length = read.infer_query_length()
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if (barcode in validated_barcodes) and (query_length >= int(min_query_length)):
            count_dic[barcode][tag][umi] += 1

    # write dic to pandas df
    rows = []
    for barcode in count_dic:
        for tag in count_dic[barcode]:
            for umi in count_dic[barcode][tag]:
                rows.append([barcode, tag, umi, count_dic[barcode][tag][umi]])
    if len(rows) == 0:
        sum_virus.logger.warning("No cell virus UMI found!")

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
        int(args.min_query_length))


def get_opts_count_capture_virus(parser, sub_program):
    parser.add_argument("--min_query_length", help='minimum query length', default=35)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='matched rna_virus directory', required=True)
        parser.add_argument('--virus_bam', required=True)

