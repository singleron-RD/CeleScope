import numpy as np
import pandas as pd

from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.vdj.__init__ import CHAINS


# UMI_min auto = CELL_CALLING_RANK UMI / 10
CELL_CALLING_RANK = 20
# mixcr sequence header
SEQUENCES_HEADER = ["aaSeqCDR3", "nSeqCDR3"]


class Count_vdj(Step):
    """
    ## Features
    - Cell-calling based on barcode-UMI rank.    
    - Summarize clonetypes infomation.

    ## Output
    - `{sample}_cell_confident.tsv` The clone type of VDJ cell barcode, each chain occupies one line.

    - `{sample}_cell_confident_count.tsv` The clone type of VDJ cell barcode, each cell occupies one line.

    - `{sample}_clonetypes.tsv` The count and percentage of each clonetypes of VDJ cell barcode.

    - `{sample}_match_clonetypes.tsv` When summarize clonetypes, only consider barcodes in the match scRNA-Seq library. 
    This file will only be produced when the `match_dir` parameter is provided.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # set
        self.chains = CHAINS[args.type]
        self.cols = []
        for chain in self.chains:
            for seq in SEQUENCES_HEADER:
                self.cols.append("_".join([seq, chain]))

        self.match_bool = False
        if utils.check_arg_not_none(self.args, 'match_dir'):
            self.match_cell_barcodes, _match_cell_number = utils.get_barcode_from_match_dir(args.match_dir)
            self.match_bool = True
        elif utils.check_arg_not_none(self.args, 'matrix_dir'):
            self.match_cell_barcodes = utils.get_barcode_from_matrix_dir(args.matrix_dir)
            self.match_bool = True
        if self.match_bool:
            self.match_cell_barcodes = set(self.match_cell_barcodes)

        # out files
        self.cell_confident_file = f"{self.out_prefix}_cell_confident.tsv"
        self.cell_confident_count_file = f"{self.out_prefix}_cell_confident_count.tsv"
        self.clonetypes_file = f"{self.out_prefix}_clonetypes.tsv"
        self.match_clonetypes_file = f"{self.out_prefix}_match_clonetypes.tsv"

        # add args data
        self.add_data(iUMI=args.iUMI)

    @utils.add_log
    def cell_calling(self, df_UMI_count_filter):
        df_UMI_sum = df_UMI_count_filter.groupby(
            ['barcode'], as_index=False).agg({"UMI": "sum"})
        if (self.args.UMI_min == "auto"):
            df_UMI_sum_sorted = df_UMI_sum.sort_values(
                ["UMI"], ascending=False)
            rank_UMI = df_UMI_sum_sorted.iloc[CELL_CALLING_RANK, :]["UMI"]
            UMI_min = int(rank_UMI / 10)
        else:
            UMI_min = int(self.args.UMI_min)
        df_UMI_cell = df_UMI_sum[df_UMI_sum.UMI >= UMI_min]
        df_UMI_sum["mark"] = df_UMI_sum["UMI"].apply(
            lambda x: "CB" if (x >= UMI_min) else "UB")

        df = df_UMI_sum.sort_values('UMI', ascending=False)
        self.add_data(CB_num=df[df['mark'] == 'CB'].shape[0])
        self.add_data(Cells=list(df.loc[df['mark'] == 'CB', 'UMI']))
        self.add_data(UB_num=df[df['mark'] == 'UB'].shape[0])
        self.add_data(Background=list(df.loc[df['mark'] == 'UB', 'UMI']))

        cell_barcodes = set(df_UMI_cell.barcode)
        total_cell_number = len(cell_barcodes)
        self.add_metric(
            name="Estimated Number of Cells",
            value=total_cell_number,
            help_info="number of barcodes considered as cell-associated"
        )

        df_cell = df_UMI_count_filter[df_UMI_count_filter.barcode.isin(
            cell_barcodes)]
        return df_cell, cell_barcodes

    @utils.add_log
    def get_df_confident(self, df_cell):
        """
        1. UMI > iUMI
        2. in chain
        """
        df_iUMI = df_cell[df_cell.UMI >= self.args.iUMI]
        df_confident = df_iUMI[df_iUMI["chain"].isin(self.chains)]
        df_confident = df_confident.sort_values("UMI", ascending=False)
        df_confident = df_confident.groupby(
            ["barcode", "chain"], as_index=False).head(1)
        return df_confident

    def get_df_valid_count(self, df_confident):
        df_valid_count = df_confident.set_index(["barcode", "chain"])
        df_valid_count = df_valid_count.unstack()
        df_valid_count.columns = ['_'.join(col) for col in df_valid_count]
        df_valid_count = df_valid_count.reset_index()
        df_valid_count.fillna(inplace=True, value="NA")
        return df_valid_count

    def get_clonetypes_and_write(self, df_valid_count, cell_barcodes):
        """
        Returns
        - df_clonetypes
        - df_match_clonetypes
        """

        total_cell_number = len(cell_barcodes)
        df_clonetypes = df_valid_count.copy()
        df_match_clonetypes = None

        df_clonetypes = df_clonetypes.groupby(self.cols, as_index=False).agg({
            "barcode": "count"})
        # put na last
        df_clonetypes.replace('NA', np.nan, inplace=True)
        df_clonetypes.sort_values(
            ["barcode"] + self.cols, ascending=False, na_position='last', inplace=True)
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
        order = ["clonetype_ID"] + self.cols + ["barcode", "percent"]
        df_clonetypes = df_clonetypes[order]
        df_clonetypes.rename(
            columns={"barcode": "barcode_count"}, inplace=True)
        # out clonetypes
        df_clonetypes.to_csv(self.clonetypes_file, sep="\t", index=False)

        if self.args.type == "TCR":

            UMI_col_dic = {"TRA": "UMI_TRA", "TRB": "UMI_TRB"}
            for chain in UMI_col_dic:
                UMI_col_name = UMI_col_dic[chain]
                if UMI_col_name in df_valid_count.columns:
                    df_valid_count[UMI_col_name].replace(
                        "NA", 0, inplace=True)
                    Median_chain_UMIs_per_Cell = np.median(
                        df_valid_count[UMI_col_name])
                else:
                    Median_chain_UMIs_per_Cell = 0
                self.add_metric(
                    name=f"Median {chain} UMIs per Cell",
                    value=Median_chain_UMIs_per_Cell,
                    help_info=f"median number of UMI mapped to {chain} per cell"
                )

            df_TRA_TRB = df_valid_count[
                (df_valid_count.aaSeqCDR3_TRA != "NA") &
                (df_valid_count.aaSeqCDR3_TRB != "NA")
            ]
            cell_with_confident_TRA_and_TRB = df_TRA_TRB.shape[0]
            self.add_metric(
                name="Cell with TRA and TRB",
                value=cell_with_confident_TRA_and_TRB,
                total=total_cell_number,
                help_info=f"cells with as least {self.args.iUMI} UMI mapped to each chain"
            )

            if self.match_bool:
                cell_with_match_barcode = self.match_cell_barcodes.intersection(
                    cell_barcodes)
                cell_with_match_barcode_number = len(cell_with_match_barcode)

                df_match = df_valid_count[df_valid_count.barcode.isin(
                    self.match_cell_barcodes)]

                df_match_TRA_TRB = df_match[
                    (df_match.aaSeqCDR3_TRA != "NA") &
                    (df_match.aaSeqCDR3_TRB != "NA")
                ]
                match_cell_with_TRA_and_TRB = df_match_TRA_TRB.shape[0]
                self.add_metric(
                    name="Cell with Barcode Match",
                    value=cell_with_match_barcode_number,
                    total=total_cell_number,
                    help_info="cells with barcode matched with scRNA-seq library"
                )
                self.add_metric(
                    name="Cell with Barcode Match, TRA and TRB",
                    value=match_cell_with_TRA_and_TRB,
                    total=cell_with_match_barcode_number,
                    help_info=f"cell with matched barcode and with as least {self.args.iUMI} UMI mapped to each chain. \
When calculating the percentage, the denominator is `Cell with Barcode Match`"
                )

        # BCR
        elif self.args.type == "BCR":

            UMI_col_dic = {"IGH": "UMI_IGH",
                           "IGL": "UMI_IGL", "IGK": "UMI_IGK"}
            for chain in UMI_col_dic:
                UMI_col_name = UMI_col_dic[chain]
                if UMI_col_name in df_valid_count.columns:
                    df_valid_count[UMI_col_name].replace(
                        "NA", 0, inplace=True)
                    df_valid_count_over_zero = df_valid_count[
                        df_valid_count[UMI_col_name] > 0
                    ]
                    Median_chain_UMIs_per_Cell = np.median(
                        df_valid_count_over_zero[UMI_col_name])
                else:
                    Median_chain_UMIs_per_Cell = 0
                self.add_metric(
                    name=f"Median {chain} UMIs per Cell",
                    value=Median_chain_UMIs_per_Cell,
                    help_info="median number of UMI mapped to each chain per cell"
                )

            df_heavy_and_light = df_valid_count[
                (df_valid_count.aaSeqCDR3_IGH != "NA") &
                (
                    (df_valid_count.aaSeqCDR3_IGL != "NA") |
                    (df_valid_count.aaSeqCDR3_IGK != "NA")
                )
            ]
            Cell_with_Heavy_and_Light_Chain = df_heavy_and_light.shape[0]
            self.add_metric(
                name="Cell with Heavy and Light Chain",
                value=Cell_with_Heavy_and_Light_Chain,
                total=total_cell_number,
                help_info=f"cells with as least {self.args.iUMI} UMI mapped to each chain"
            )

            if self.match_bool:
                cell_with_match_barcode = self.match_cell_barcodes.intersection(
                    cell_barcodes)
                cell_with_match_barcode_number = len(cell_with_match_barcode)

                df_match = df_valid_count[df_valid_count.barcode.isin(
                    self.match_cell_barcodes)]

                # median match UMI
                df_match_heavy_light = df_match[
                    (df_match.aaSeqCDR3_IGH != "NA") &
                    (
                        (df_match.aaSeqCDR3_IGL != "NA") |
                        (df_match.aaSeqCDR3_IGK != "NA")
                    )
                ]
                match_cell_with_heavy_and_light = df_match_heavy_light.shape[0]
                self.add_metric(
                    name="Cell with Barcode Match",
                    value=cell_with_match_barcode_number,
                    total=total_cell_number,
                    help_info="cells with barcode matched with scRNA-seq library"
                )
                self.add_metric(
                    name="Cell with Barcode Match, Heavy and Light Chain",
                    value=match_cell_with_heavy_and_light,
                    total=total_cell_number,
                    help_info=f"cell with matched barcode and with as least {self.args.iUMI} UMI mapped to each chain"
                )

        if self.match_bool:
            """
            df_match_clonetypes
            """
            df_match_clonetypes = df_match.groupby(self.cols, as_index=False).agg({
                "barcode": "count"})
            total_match_CDR3_barcode_number = sum(
                df_match_clonetypes.barcode)
            df_match_clonetypes["percent"] = df_match_clonetypes.barcode / \
                total_match_CDR3_barcode_number * 100
            df_match_clonetypes["percent"] = df_match_clonetypes["percent"].apply(
                lambda x: round(x, 2)
            )
            df_match_clonetypes.rename(
                columns={"barcode": "barcode_count"}, inplace=True)
            df_match_clonetypes = df_match_clonetypes.merge(
                df_clonetypes, on=self.cols, how='left', suffixes=('', '_y'))
            # order and drop duplicated cols
            order = ["clonetype_ID"] + self.cols + ["barcode_count", "percent"]
            df_match_clonetypes = df_match_clonetypes[order]
            df_match_clonetypes.sort_values(["barcode_count", "clonetype_ID"], ascending=[
                                            False, True], inplace=True)
            df_match_clonetypes.to_csv(
                self.match_clonetypes_file, sep="\t", index=False)
        return df_clonetypes, df_match_clonetypes

    def write_cell_confident_count(self, df_valid_count, df_clonetypes, df_confident):
        df_mergeID = pd.merge(df_valid_count,
                              df_clonetypes, how="left", on=self.cols)
        df_mergeID.sort_values(["clonetype_ID", "barcode"], inplace=True)
        # output df_valid_count
        df_mergeID.to_csv(self.cell_confident_count_file,
                          sep="\t", index=False)
        df_mergeID = df_mergeID[["barcode", "clonetype_ID"]]
        df_cell_confident_with_ID = pd.merge(
            df_confident, df_mergeID, how="left", on="barcode")
        df_cell_confident_with_ID.sort_values(
            ["clonetype_ID", "barcode", "chain"], inplace=True)
        # output df_cell_confident
        df_cell_confident_with_ID.to_csv(
            self.cell_confident_file, sep="\t", index=False)

    def write_clonetypes_table_to_data(self, df_clonetypes, df_match_clonetypes):
        # cloneytpes table
        def format_table(df_clonetypes):
            df_table = df_clonetypes.copy()
            df_table["percent"] = df_table["percent"].apply(
                lambda x: str(x) + "%")
            seqs = ["aaSeqCDR3"]
            cols = []
            for chain in self.chains:
                for seq in seqs:
                    cols.append("_".join([seq, chain]))
            df_table_cols = ["clonetype_ID"] + \
                cols + ["barcode_count", "percent"]
            df_table = df_table[df_table_cols]
            table_header = ["Clonetype_ID"] + cols + ["Frequency", "Percent"]
            return df_table, table_header

        df_table, _table_header = format_table(df_clonetypes)
        title = 'Clonetypes'
        if self.match_bool:
            df_table, _table_header = format_table(df_match_clonetypes)
            title = 'Match Clonetypes'

        table_dict = self.get_table_dict(
            title=title,
            table_id='clonetypes',
            df_table=df_table
        )
        self.add_data(table_dict=table_dict)

    def run(self):
        df_UMI_count_filter = pd.read_csv(
            self.args.UMI_count_filter_file, sep='\t')
        df_cell, cell_barcodes = self.cell_calling(df_UMI_count_filter)
        df_confident = self.get_df_confident(df_cell)
        df_valid_count = self.get_df_valid_count(df_confident)
        df_clonetypes, df_match_clonetypes = self.get_clonetypes_and_write(
            df_valid_count, cell_barcodes)
        self.write_cell_confident_count(
            df_valid_count, df_clonetypes, df_confident)
        self.write_clonetypes_table_to_data(df_clonetypes, df_match_clonetypes)


def count_vdj(args):
    # TODO
    # add TCR or BCR prefix to distinguish them in html report summary; should improve
    with Count_vdj(args, display_title="Count") as runner:
        runner.run()


def get_opts_count_vdj(parser, sub_program):
    parser.add_argument("--type", help=HELP_DICT['type'], required=True)
    parser.add_argument('--UMI_min',help=HELP_DICT['UMI_min'],default="auto")
    parser.add_argument(
        '--iUMI',
        help="""Default `1`. Minimum number of UMI of identical receptor type and CDR3. 
For each (barcode, chain) combination, only UMI>=iUMI is considered valid.""",
        type=int,
        default=1
    )
    if sub_program:
        parser.add_argument("--UMI_count_filter_file",help=HELP_DICT['UMI_count_filter_file'], required=True)
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--matrix_dir", help=HELP_DICT['matrix_dir'])
        parser = s_common(parser)
