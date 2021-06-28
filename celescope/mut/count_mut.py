import os

import editdistance
import pandas as pd
import pysam

from celescope.tools.utils import add_log, parse_match_dir

parentDir = os.path.dirname(__file__)


def read_mut(mut_file):
    mut_dic = {}
    df = pd.read_csv(mut_file, sep="\t")

    def convert_row(row):
        gene = row["gene"]
        mut_type = row["type"]
        seq = row["seq"]
        ref_position = row["ref_position"]
        if gene not in mut_dic:
            mut_dic[gene] = {}
        mut_dic[gene][mut_type] = {}
        mut_dic[gene][mut_type]["pos"] = ref_position
        mut_dic[gene][mut_type]["seq"] = seq
    df.apply(func=convert_row, axis=1)
    return mut_dic


def is_fusion(pos, read_start, read_length, flanking_base):
    test_start = (pos - flanking_base) >= read_start
    test_end = (pos + flanking_base) <= (read_start + read_length)
    return (test_start and test_end)


@add_log
def count_mut(args):

    outdir = args.outdir
    sample = args.sample
    bam = args.bam
    shift_base = int(args.shift_base)
    mut_file = args.mut_file
    match_dir = args.match_dir

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    mut_dic = read_mut(mut_file)
    out_prefix = outdir + "/" + sample

    # tsne
    match_dict = parse_match_dir(match_dir)
    df_tsne = pd.read_csv(match_dict['tsne_coord'], sep="\t", index_col=0)

    # out
    out_read_file = out_prefix + "_mut_read.tsv"
    out_read_count_file = out_prefix + "_mut_read_count.tsv"
    out_umi_count_file = out_prefix + "_mut_UMI_count.tsv"
    out_insertion_barcode_count_file = out_prefix + "_mut_insertion_barcode_count.tsv"
    out_tsne_file = out_prefix + "_mut_tsne.tsv"

    # process bam
    samfile = pysam.AlignmentFile(bam, "rb")
    # header = samfile.header
    # new_bam = pysam.AlignmentFile(out_prefix+"_mut.bam", "wb", header=header)
    rows = []
    for read in samfile:
        tag = read.reference_name
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        seq = read.query_sequence
        cigar = read.cigar
        ref_pos = read.reference_start
        read_pos = read.query_alignment_start
        for (cigarType, cigarLength) in cigar:
            if cigarType == 0:  # match
                ref_pos += cigarLength
                read_pos += cigarLength
            elif cigarType == 1:  # insertion
                insert_seq = seq[read_pos:read_pos+cigarLength]
                insert_seq_length = len(insert_seq)
                rows.append({
                    "barcode": barcode,
                    "UMI": umi,
                    "gene": tag,
                    "type": "insertion",
                    "seq": insert_seq,
                    "ref_pos": ref_pos,
                    "read_pos": read_pos,
                    "seq_length": insert_seq_length
                })
                read_pos += cigarLength
            elif cigarType == 2:  # deletion
                ref_pos += cigarLength
    df = pd.DataFrame(rows)
    df = df[["barcode", "UMI", "gene", "type",
             "ref_pos", "read_pos", "seq", "seq_length"]]

    def is_valid(row):
        gene = str(row["gene"])
        ref_pos = int(row["ref_pos"])
        mut_type = str(row["type"])
        int(row["seq_length"])
        seq = str(row["seq"])
        if gene not in mut_dic:
            return False
        if mut_type not in mut_dic[gene].keys():
            return False
        if abs(ref_pos-mut_dic[gene][mut_type]["pos"]) > shift_base:
            return False
        if editdistance.eval(seq, mut_dic[gene][mut_type]["seq"]) > shift_base:
            return False
        return True

    df["is_valid"] = df.apply(func=is_valid, axis=1)
    df.to_csv(out_read_file, sep="\t")

    df_valid = df[df["is_valid"]]
    df_read_count = df_valid.groupby(
        ["gene", "type", "barcode", "UMI"]).agg({"UMI": "count"})
    df_read_count.columns = ["read_count"]
    df_read_count.reset_index(inplace=True)
    df_read_count.to_csv(out_read_count_file, sep="\t")

    df_UMI_count = df_read_count.groupby(
        ["gene", "type", "barcode"]).agg({"UMI": "count"})
    df_UMI_count.columns = ["UMI_count"]
    df_UMI_count.to_csv(out_umi_count_file, sep="\t")

    df_temp = df_UMI_count.reset_index()
    if df_temp.shape[0] == 0:
        count_mut.logger.warning('NO VALID INSERTION FOUND!')
    else:
        df_insertion = df_temp[df_temp["type"] == "insertion"]
        df_insertion_barcode_count = df_insertion.pivot(
            index="barcode", columns="gene", values="UMI_count")
        df_insertion_barcode_count.to_csv(
            out_insertion_barcode_count_file, sep="\t")

        df_tsne_mut = pd.merge(df_tsne, df_insertion_barcode_count,
                               right_index=True, left_index=True, how="left")
        df_tsne_mut.fillna(0, inplace=True)
        df_tsne_mut.to_csv(out_tsne_file, sep="\t")

        # plot
        app = parentDir + "/plot.R"
        cmd = f"Rscript {app} --tsne_mut {out_tsne_file} --outdir {outdir}"
        count_mut.logger.info(cmd)
        os.system(cmd)
        count_mut.logger.info("plot done.")


def get_opts_count_mut(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--bam", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument("--mut_file", help="mutation file", required=True)
    parser.add_argument(
        "--match_dir", help="match scRNA-Seq dir", required=True)
    parser.add_argument("--shift_base", default=2)
