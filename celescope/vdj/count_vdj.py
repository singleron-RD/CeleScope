import glob
from celescope.vdj.__init__ import CHAINS
from celescope.tools.report import reporter
from celescope.tools.utils import format_number, log, read_barcode_file
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
mpl.use('Agg')
from matplotlib import pyplot as plt


def report_prepare(df, outdir):

    json_file = outdir + '/.data.json'
    if not os.path.exists(json_file):
        data = {}
    else:
        fh = open(json_file)
        data = json.load(fh)
        fh.close()

    df = df.sort_values('UMI', ascending=False)
    data['CB_num'] = df[df['mark'] == 'CB'].shape[0]
    data['Cells'] = list(df.loc[df['mark'] == 'CB', 'UMI'])
    data['UB_num'] = df[df['mark'] == 'UB'].shape[0]
    data['Background'] = list(df.loc[df['mark'] == 'UB', 'UMI'])

    with open(json_file, 'w') as fh:
        json.dump(data, fh)


@log
def count_vdj(args):

    sample = args.sample
    match_dir = args.match_dir
    UMI_min = args.UMI_min
    outdir = args.outdir
    UMI_count_filter1_file = args.UMI_count_filter1_file
    type = args.type
    debug = args.debug
    iUMI = int(args.iUMI)
    chains = CHAINS[type]

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % outdir)

    # out file
    cell_confident_file = f"{outdir}/{sample}_cell_confident.tsv"
    cell_confident_count_file = f"{outdir}/{sample}_cell_confident_count.tsv"
    clonetypes_file = f"{outdir}/{sample}_clonetypes.tsv"
    match_clonetypes_file = f"{outdir}/{sample}_match_clonetypes.tsv"
    top10_clonetypes_file = f"{outdir}/{sample}_top10_clonetypes.tsv"
    match_top10_clonetypes_file = f"{outdir}/{sample}_match_top10_clonetypes.tsv"

    # read file
    df_UMI_count_filter1 = pd.read_csv(UMI_count_filter1_file, sep='\t')
    if (not match_dir) or (match_dir == "None"):
        match_bool = False
    else:
        match_bool = True
    if match_bool:
        match_cell_barcodes, match_cell_number = read_barcode_file(match_dir)

    cell_summary_row_list = []

    # cell calling：cell calling: keep UMIs >= UMI_min
    df_UMI_sum = df_UMI_count_filter1.groupby(
        ['barcode'], as_index=False).agg({"UMI": "sum"})
    if (UMI_min == "auto"):
        rank = 20
        df_UMI_sum_sorted = df_UMI_sum.sort_values(["UMI"], ascending=False)
        rank_UMI = df_UMI_sum_sorted.iloc[rank, :]["UMI"]
        UMI_min = int(rank_UMI / 10)
    else:
        UMI_min = int(UMI_min)
    df_UMI_cell = df_UMI_sum[df_UMI_sum.UMI >= UMI_min]
    df_UMI_sum["mark"] = df_UMI_sum["UMI"].apply(
        lambda x: "CB" if (x >= UMI_min) else "UB")
    report_prepare(df_UMI_sum, outdir + "/../")

    cell_barcodes = set(df_UMI_cell.barcode)
    cell_number = len(cell_barcodes)
    cell_summary_row_list.append({
        "item": "Estimated Number of Cells",
        "count": cell_number,
        "total_count": cell_number,
    })

    # df_UMI_count_filter1 in cell
    df_cell = df_UMI_count_filter1[df_UMI_count_filter1.barcode.isin(
        cell_barcodes)]
    # filter2: cell wtih UMI >= iUMI of identical receptor type and CDR3
    # combinations.
    df_cell_UMI_count_filter2 = df_cell[df_cell.UMI >= iUMI]

    # cell confident
    df_cell_confident = df_cell_UMI_count_filter2[df_cell_UMI_count_filter2["chain"].isin(
        chains)]
    df_cell_confident = df_cell_confident.sort_values("UMI", ascending=False)
    df_cell_confident = df_cell_confident.groupby(
        ["barcode", "chain"], as_index=False).head(1)

    # count
    df_cell_confident_count = df_cell_confident.set_index(["barcode", "chain"])
    df_cell_confident_count = df_cell_confident_count.unstack()
    df_cell_confident_count.columns = [
        '_'.join(col) for col in df_cell_confident_count]
    df_cell_confident_count = df_cell_confident_count.reset_index()
    df_cell_confident_count.fillna(inplace=True, value="NA")

    # clonetypes
    seqs = ["aaSeqCDR3", "nSeqCDR3"]
    cols = []
    for chain in chains:
        for seq in seqs:
            cols.append("_".join([seq, chain]))

    for col in cols:
        if not (col in list(df_cell_confident_count.columns)):
            df_cell_confident_count[col] = "NA"

    df_clonetypes = df_cell_confident_count.copy()

    df_clonetypes = df_clonetypes.groupby(cols, as_index=False).agg({
        "barcode": "count"})
    # put na last
    df_clonetypes.replace('NA', np.nan, inplace=True)
    df_clonetypes.sort_values(["barcode"] + cols, ascending=False, na_position='last', inplace=True)
    df_clonetypes.replace(np.nan, 'NA', inplace=True)

    total_CDR3_barcode_number = sum(df_clonetypes.barcode)
    df_clonetypes["percent"] = df_clonetypes.barcode / \
        total_CDR3_barcode_number * 100
    df_clonetypes["percent"] = df_clonetypes["percent"].apply(
        lambda x: round(x, 2))

    # add clonetype ID
    df_clonetypes = df_clonetypes.reset_index()
    df_clonetypes["clonetype_ID"] = pd.Series(df_clonetypes.index) + 1
    df_clonetypes.drop(columns=["index"], inplace=True)

    # order
    order = ["clonetype_ID"] + cols + ["barcode", "percent"]
    df_clonetypes = df_clonetypes[order]
    df_clonetypes.rename(columns={"barcode": "barcode_count"}, inplace=True)
    # out clonetypes
    df_clonetypes.to_csv(clonetypes_file, sep="\t", index=False)

    if type == "TCR":

        UMI_col_dic = {"TRA": "UMI_TRA", "TRB": "UMI_TRB"}
        for chain in UMI_col_dic:
            UMI_col_name = UMI_col_dic[chain]
            if UMI_col_name in df_cell_confident_count.columns:
                df_cell_confident_count[UMI_col_name].replace(
                    "NA", 0, inplace=True)
                Median_chain_UMIs_per_Cell = np.median(
                    df_cell_confident_count[UMI_col_name])
            else:
                Median_chain_UMIs_per_Cell = 0
            cell_summary_row_list.append({
                "item": "Median {chain} UMIs per Cell".format(chain=chain),
                "count": Median_chain_UMIs_per_Cell,
                "total_count": np.nan
            })

        df_TRA_TRB = df_cell_confident_count[
            (df_cell_confident_count.aaSeqCDR3_TRA != "NA") &
            (df_cell_confident_count.aaSeqCDR3_TRB != "NA")
        ]
        cell_with_confident_TRA_and_TRB = df_TRA_TRB.shape[0]
        cell_summary_row_list.append({
            "item": "Cell with TRA and TRB",
            "count": cell_with_confident_TRA_and_TRB,
            "total_count": cell_number,
        })

        """
        df cell barcode filter
        intersect cell_barcodes from scRNA-Seq with barcode from TCR seq
        """
        if match_bool:
            cell_with_match_barcode = match_cell_barcodes.intersection(
                cell_barcodes)
            cell_with_match_barcode_number = len(cell_with_match_barcode)

            df_match = df_cell_confident_count[df_cell_confident_count.barcode.isin(
                match_cell_barcodes)]

            df_match_TRA_TRB = df_match[
                (df_match.aaSeqCDR3_TRA != "NA") &
                (df_match.aaSeqCDR3_TRB != "NA")
            ]
            match_cell_with_TRA_and_TRB = df_match_TRA_TRB.shape[0]

            cell_summary_row_list.append({
                "item": "Cell with Barcode Match",
                "count": cell_with_match_barcode_number,
                "total_count": cell_number,
            })
            cell_summary_row_list.append({
                "item": "Cell with Barcode Match, TRA and TRB",
                "count": match_cell_with_TRA_and_TRB,
                "total_count": cell_number,
            })

    # BCR
    elif type == "BCR":

        UMI_col_dic = {"IGH": "UMI_IGH", "IGL": "UMI_IGL", "IGK": "UMI_IGK"}
        for chain in UMI_col_dic:
            UMI_col_name = UMI_col_dic[chain]
            if UMI_col_name in df_cell_confident_count.columns:
                df_cell_confident_count[UMI_col_name].replace(
                    "NA", 0, inplace=True)
                df_cell_confident_count_over_zero = df_cell_confident_count[
                    df_cell_confident_count[UMI_col_name] > 0
                ]
                Median_chain_UMIs_per_Cell = np.median(
                    df_cell_confident_count_over_zero[UMI_col_name])
            else:
                Median_chain_UMIs_per_Cell = 0
            cell_summary_row_list.append({
                "item": "Median {chain} UMIs per Cell".format(chain=chain),
                "count": Median_chain_UMIs_per_Cell,
                "total_count": np.nan})

        df_heavy_and_light = df_cell_confident_count[
            (df_cell_confident_count.aaSeqCDR3_IGH != "NA") &
            (
                (df_cell_confident_count.aaSeqCDR3_IGL != "NA") |
                (df_cell_confident_count.aaSeqCDR3_IGK != "NA")
            )
        ]
        Cell_with_Heavy_and_Light_Chain = df_heavy_and_light.shape[0]
        cell_summary_row_list.append({
            "item": "Cell with Heavy and Light Chain",
            "count": Cell_with_Heavy_and_Light_Chain,
            "total_count": cell_number
        })

        """
        df cell barcode filter
        intersect cell_barcodes from normal scRNA-Seq with barcode from BCR seq
        """
        if match_bool:
            cell_with_match_barcode = match_cell_barcodes.intersection(
                cell_barcodes)
            cell_with_match_barcode_number = len(cell_with_match_barcode)

            df_match = df_cell_confident_count[df_cell_confident_count.barcode.isin(
                match_cell_barcodes)]

            # median match UMI
            df_match_heavy_light = df_match[
                (df_match.aaSeqCDR3_IGH != "NA") &
                (
                    (df_match.aaSeqCDR3_IGL != "NA") |
                    (df_match.aaSeqCDR3_IGK != "NA")
                )
            ]
            match_cell_with_heavy_and_light = df_match_heavy_light.shape[0]

            cell_summary_row_list.append({
                "item": "Cell with Barcode Match ",
                "count": cell_with_match_barcode_number,
                "total_count": cell_number
            })
            cell_summary_row_list.append({
                "item": "Cell with Barcode Match, Heavy and Light Chain",
                "count": match_cell_with_heavy_and_light,
                "total_count": cell_number
            })

    if match_bool:
        """
        df_match_clonetypes
        """
        df_match_clonetypes = df_match.groupby(cols, as_index=False).agg({
            "barcode": "count"})
        total_match_CDR3_barcode_number = sum(
            df_match_clonetypes.barcode)
        df_match_clonetypes["percent"] = df_match_clonetypes.barcode / \
            total_match_CDR3_barcode_number * 100
        df_match_clonetypes["percent"] = df_match_clonetypes["percent"].apply(
            lambda x: round(x, 2)
        )
        df_match_clonetypes.rename(columns={"barcode": "barcode_count"}, inplace=True)
        df_match_clonetypes = df_match_clonetypes.merge(
            df_clonetypes, on=cols, how='left', suffixes=('', '_y'))
        # order and drop duplicated cols
        order = ["clonetype_ID"] + cols + ["barcode_count", "percent"]
        df_match_clonetypes = df_match_clonetypes[order]
        df_match_clonetypes.sort_values(["barcode_count", "clonetype_ID"], ascending=[False,True], inplace=True)
        df_match_clonetypes.to_csv(
            match_clonetypes_file, sep="\t", index=False)

    df_mergeID = pd.merge(df_cell_confident_count,
                          df_clonetypes, how="left", on=cols)
    df_mergeID.sort_values(["clonetype_ID", "barcode"], inplace=True)
    # output df_cell_confident_count
    df_mergeID.to_csv(cell_confident_count_file, sep="\t", index=False)
    df_mergeID = df_mergeID[["barcode", "clonetype_ID"]]
    df_cell_confident_with_ID = pd.merge(
        df_cell_confident, df_mergeID, how="left", on="barcode")
    df_cell_confident_with_ID.sort_values(
        ["clonetype_ID", "barcode", "chain"], inplace=True)
    # output df_cell_confident
    df_cell_confident_with_ID.to_csv(
        cell_confident_file, sep="\t", index=False)

    # summary file
    cell_summary = pd.DataFrame(cell_summary_row_list, columns=[
                                "item", "count", "total_count"])
    cell_summary["count"] = cell_summary["count"].apply(int)
    cell_summary["percent"] = cell_summary["count"] / \
        (cell_summary.total_count.astype("float")) * 100
    cell_summary["percent"] = cell_summary["percent"].apply(
        lambda x: round(x, 2))
    cell_summary["count"] = cell_summary["count"].apply(format_number)

    def percent_str_func(row):
        need_percent = bool(
            re.search("Cell with", row["item"], flags=re.IGNORECASE))
        if need_percent:
            return "(" + str(row["percent"]) + "%)"
        else:
            return ""
    cell_summary["percent_str"] = cell_summary.apply(
        lambda row: percent_str_func(row), axis=1)

    # stat file
    def gen_stat(summary, stat_file):
        stat = summary
        stat["new_count"] = stat["count"].astype(str) + stat["percent_str"]
        stat = stat.loc[:, ["item", "new_count"]]
        stat.to_csv(stat_file, sep=":", header=None, index=False)

    cell_stat_file = "{}/stat.txt".format(outdir)
    gen_stat(cell_summary, cell_stat_file)
    name = type + '_count_vdj'
    t = reporter(
        name=name,
        sample=args.sample,
        stat_file=cell_stat_file,
        outdir=outdir + '/..',
        assay=args.assay,
        parameters={"iUMI": iUMI},
    )
    t.get_report()

    # cloneytpes table
    def format_table(df_clonetypes, top10_clonetypes_file):
        top10_clonetypes_df = df_clonetypes.head(10)
        top10_clonetypes_df = top10_clonetypes_df.reset_index(drop=True)
        top10_clonetypes_df.index = top10_clonetypes_df.index + 1
        top10_clonetypes_df["percent"] = top10_clonetypes_df["percent"].apply(
            lambda x: str(x) + "%")
        seqs = ["aaSeqCDR3"]
        cols = []
        for chain in chains:
            for seq in seqs:
                cols.append("_".join([seq, chain]))
        top10_cols = ["clonetype_ID"] + cols + ["barcode_count", "percent"]
        top10_clonetypes_df = top10_clonetypes_df[top10_cols]
        top10_clonetypes_df.to_csv(top10_clonetypes_file, sep="\t", index=False)
        table_header = ["Clonetype_ID"] + cols + ["Frequency", "Percent"]
        return table_header

    table_header = format_table(df_clonetypes, top10_clonetypes_file)
    use_top10_clonetypes_file = top10_clonetypes_file
    section_header = 'Top10 clonetypes'
    if match_bool:
        format_table(df_match_clonetypes, match_top10_clonetypes_file)
        use_top10_clonetypes_file = match_top10_clonetypes_file
        section_header = 'Match Top10 clonetypes'
    
    t = reporter(
        name="clonetypes",
        sample=args.sample,
        table_file=use_top10_clonetypes_file,
        table_header=table_header,
        outdir=outdir + '/..',
        assay=args.assay,
        parameters={'section_header': section_header},
    )
    t.get_report()

    # other_metrics_file
    """
    if len(other_metrics_row_list) != 0:
        other_metrics = pd.DataFrame(other_metrics_row_list,columns=["item","count"])
        other_metrics.to_csv(other_metrics_file,sep=":",header=None,index=False)
    """


def get_opts_count_vdj(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--UMI_count_filter1_file", required=True)
        parser.add_argument("--assay", required=True)
    parser.add_argument("--match_dir", default=None)
    parser.add_argument("--type", required=True)
    parser.add_argument('--UMI_min', dest='UMI_min',
                        help='minimum UMI number to filter', default="auto")
    parser.add_argument('--debug', dest='debug', default=False)
    parser.add_argument(
        '--iUMI', help='minimum number of UMI of identical receptor type and CDR3', default=1)
