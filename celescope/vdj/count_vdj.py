import numpy as np
import pandas as pd

from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.vdj.__init__ import CHAINS, PAIRS
from celescope.tools.capture.threshold import Auto
from celescope.tools.emptydrop_cr import get_plot_elements

# mixcr sequence header
SEQUENCES_HEADER = ["aaSeqCDR3", "nSeqCDR3"]
# Auto filter percentitle
PERCENTILE = {
    'TCR': 85,
    'BCR': 90,
}


def target_cell_calling(df_UMI_sum, expected_target_cell_num=3000, target_barcodes=None, weight=3, coef=10, 
    percentile=99, UMI_min=3, umi_col='UMI'):
    """
    Args:
        df_UMI_sum: A dataframe with columns barcode and UMI. The barcode are match cell barcodes.
    
    Returns:
        umi_threhold: int
        target_cell_barcodes: list

    >>> df_UMI_sum = pd.DataFrame({"barcode": ["A", "B", "C", "D", "E"], "UMI": [1, 2, 1, 30, 40]})
    >>> umi_threshold, target_cell_barcodes = target_cell_calling(df_UMI_sum, expected_target_cell_num=5, percentile=80, coef=10, target_barcodes=["A", "C"])
    >>> umi_threshold == 3
    True
    >>> target_cell_barcodes == {'A', 'C', 'D', 'E'}
    True
    """
    umi_threshold = Auto(list(df_UMI_sum[umi_col]), expected_cell_num=expected_target_cell_num, coef=coef, percentile=percentile).run()
    if umi_threshold < UMI_min:
        umi_threshold = UMI_min
    # avoid change the original dataframe
    df_temp = df_UMI_sum.copy()
    if target_barcodes:
        df_temp[umi_col] = df_temp.apply(
            lambda row:  row[umi_col] * weight if row['barcode'] in target_barcodes else row[umi_col], axis=1)
             
    target_cell_barcodes = set(df_temp.loc[df_temp[umi_col] >= umi_threshold].barcode)

    return umi_threshold, target_cell_barcodes


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
        self.pairs = PAIRS[args.type]
        self.cols = []
        for chain in self.chains:
            for seq in SEQUENCES_HEADER:
                self.cols.append("_".join([seq, chain]))

        if utils.check_arg_not_none(self.args, 'match_dir'):
            self.match_cell_barcodes, _match_cell_number = utils.get_barcode_from_match_dir(args.match_dir)
        elif utils.check_arg_not_none(self.args, 'matrix_dir'):
            self.match_cell_barcodes, _match_cell_number = utils.get_barcode_from_matrix_dir(args.matrix_dir)
        else:
            raise FileNotFoundError("--match_dir or --matrix_dir is required.")
        self.match_cell_barcodes = set(self.match_cell_barcodes)

        # input data
        df_UMI_count_filter = pd.read_csv(args.UMI_count_filter_file, sep='\t')
        self.df_UMI_sum = df_UMI_count_filter.groupby(
            ['barcode'], as_index=False).agg({"UMI": "sum"})
        self.df_match_UMI_count_filter = df_UMI_count_filter[df_UMI_count_filter['barcode'].isin(self.match_cell_barcodes)]
        self.df_match_UMI_sum = self.df_UMI_sum[self.df_UMI_sum['barcode'].isin(self.match_cell_barcodes)]

        if args.target_cell_barcode:
            self.target_barcodes = utils.read_one_col(args.target_cell_barcode)
            self.expected_target_cell_num = len(self.target_barcodes)
        else:
            self.target_barcodes = None
            self.expected_target_cell_num = args.expected_target_cell_num


        # out files
        self.UMI_sum_file = f"{self.out_prefix}_UMI_sum.tsv"
        self.cell_confident_file = f"{self.out_prefix}_cell_confident.tsv"
        self.cell_confident_count_file = f"{self.out_prefix}_cell_confident_count.tsv"
        self.clonetypes_file = f"{self.out_prefix}_clonetypes.tsv"
        self.match_clonetypes_file = f"{self.out_prefix}_match_clonetypes.tsv"
        self.filtered_annotation = f"{self.out_prefix}_filtered_contig_annotations.csv"


    @utils.add_log
    def cell_calling(self):
        """

        """
        # UMI_min filter
        percentile = PERCENTILE[self.args.type]
        umi_threshold, target_cell_barcodes = target_cell_calling(
            self.df_match_UMI_sum, 
            expected_target_cell_num=self.expected_target_cell_num, 
            target_barcodes=self.target_barcodes,
            UMI_min=self.args.UMI_min,
            percentile=percentile
        )
        df_cell = self.df_match_UMI_count_filter[self.df_match_UMI_count_filter['barcode'].isin(
            target_cell_barcodes)]
        self.df_UMI_sum['mark'] = self.df_UMI_sum['barcode'].apply(lambda x: 'CB' if x in target_cell_barcodes else 'UB')
        self.df_UMI_sum = self.df_UMI_sum.sort_values(by=['UMI'], ascending=False)
        self.df_UMI_sum.to_csv(self.UMI_sum_file, sep='\t', index=False)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.UMI_sum_file))

        self.add_metric(
            name="Estimated Number of Cells",
            value=len(target_cell_barcodes),
            help_info=f"number of barcodes considered as cell-associated. \
                UMI_threshold: {umi_threshold}. expected_cell_num: {self.expected_target_cell_num}"
        )

        return df_cell, target_cell_barcodes

    @utils.add_log
    def get_df_confident(self, df_cell):
        """
        1. UMI > iUMI
        2. in chain
        """
        if self.args.type == 'TCR':
            iUMI = self.args.TCR_iUMI
        elif self.args.type == 'BCR':
            iUMI = self.args.BCR_iUMI
        df_iUMI = df_cell[df_cell.UMI >= iUMI]
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

    def get_clonetypes_and_write(self, df_valid_count):
        """
        Returns
        - df_clonetypes
        """

        df_clonetypes = df_valid_count.copy()
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

        return df_clonetypes

    def add_metrics(self, df_valid_count, cell_barcodes):
        n_cell = len(cell_barcodes)

        for chain in self.chains:
            UMI_col_name = "UMI_" + chain
            if UMI_col_name in df_valid_count.columns:
                df_valid_count[UMI_col_name].replace("NA", 0, inplace=True)
                Median_chain_UMIs_per_Cell = np.median(df_valid_count[UMI_col_name])
            else:
                Median_chain_UMIs_per_Cell = 0

            self.add_metric(
                name=f"Median {chain} UMIs per Cell",
                value=Median_chain_UMIs_per_Cell,
                help_info=f"median number of UMI mapped to {chain} per cell"
            )

        if self.args.type == 'TCR':
            iUMI = self.args.TCR_iUMI
        elif self.args.type == 'BCR':
            iUMI = self.args.BCR_iUMI

        for pair in self.pairs:
            df_pair = df_valid_count.copy()
            for chain in pair:
                cdr3_col_name = "aaSeqCDR3_" + chain
                df_pair = df_pair[df_pair[cdr3_col_name] != "NA"]
            
            n_cell_pair = len(df_pair.barcode.unique())

            pair_str = ','.join(pair)
            self.add_metric(
                name=f"Cells with ({pair_str}) pair",
                value=n_cell_pair,
                total=n_cell,
                help_info=f"cells with as least {iUMI} UMI mapped to each chain"
            )

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

    def write_clonetypes_table_to_data(self, df_clonetypes):
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

        table_dict = self.get_table_dict(
            title=title,
            table_id='clonetypes',
            df_table=df_table
        )
        self.add_data(table_dict=table_dict)
    
    @utils.add_log
    def convert_to_annotation(self):
        in_file = pd.read_csv(self.cell_confident_file, sep='\t')
        out_file = pd.DataFrame({
                        'barcode':in_file.barcode,
                        'chain':in_file.chain,
                        'v_gene':in_file.bestVGene,
                        'd_gene':'None',
                        'j_gene':in_file.bestJGene,
                        'c_gene':'None',
                        'full_length':True,
                        'productive':True,
                        'cdr3':in_file.aaSeqCDR3,
                        'cdr3_nt':in_file.nSeqCDR3,
                        'reads':in_file.UMI,
                        'umis':in_file.UMI,
                        'raw_clonotype_id':in_file.clonetype_ID
                        })
        out_file.to_csv(self.filtered_annotation, sep=',', index=False)

    def run(self):
        df_cell, cell_barcodes = self.cell_calling()
        df_confident = self.get_df_confident(df_cell)
        df_valid_count = self.get_df_valid_count(df_confident)
        df_clonetypes= self.get_clonetypes_and_write(
            df_valid_count)
        self.add_metrics(df_valid_count, cell_barcodes)
        self.write_cell_confident_count(
            df_valid_count, df_clonetypes, df_confident)
        self.write_clonetypes_table_to_data(df_clonetypes)
        self.convert_to_annotation()


def count_vdj(args):
    with Count_vdj(args, display_title="Count") as runner:
        runner.run()


def get_opts_count_vdj(parser, sub_program):
    parser.add_argument(
        "--type", help="Required. `TCR` or `BCR`. ", required=True)
    parser.add_argument(
        '--UMI_min',
        help='minimum number of chain UMI to consider as as cell',
        type=int,
        default=3,
    )
    parser.add_argument(
        '--BCR_iUMI',
        help="""Minimum number of UMI of identical receptor type and CDR3 for BCR. 
For each (barcode, chain) combination, only UMI>=iUMI is considered valid.""",
        type=int,
        default=2,
    )
    parser.add_argument(
        '--TCR_iUMI',
        help="""Minimum number of UMI of identical receptor type and CDR3 for BCR. 
For each (barcode, chain) combination, only UMI>=iUMI is considered valid.""",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--expected_target_cell_num", 
        help="Expected T or B cell number. If `--target_cell_barcode` is provided, this argument is ignored.", 
        type=int,
        default=3000,
    )
    parser.add_argument(
        "--target_cell_barcode", 
        help="Barcode of target cells. It is a plain text file with one barcode per line. If provided, `--expected_target_cell_num` is ignored.",
        default=None,
    )
    parser.add_argument(
        "--target_weight", 
        help="UMIs of the target cells are multiplied by this factor. Only used when `--target_cell_barcode` is provided.", 
        type=float,
        default=3.0,
    )
    if sub_program:
        parser.add_argument("--UMI_count_filter_file",
                            help="Required. File from step mapping_vdj.", required=True)
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--matrix_dir", help=HELP_DICT['matrix_dir'])

        parser = s_common(parser)
