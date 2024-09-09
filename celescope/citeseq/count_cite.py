import sys

import numpy as np
import pandas as pd

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT
from celescope.tools.matrix import CountMatrix, Features
from celescope.tools.analysis_wrapper import read_tsne


TAG_COL = "tag_name"


class Count_cite(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        self.df_read_count = pd.read_csv(
            args.read_count_file, sep="\t", index_col=[0, 1]
        )

        self.match_dict = utils.parse_match_dir(args.match_dir)
        self.match_barcode = self.match_dict["match_barcode"]
        self.match_matrix_dir = self.match_dict["matrix_dir"]
        self.df_rna_tsne = pd.DataFrame()
        if "tsne_coord" in self.match_dict:
            self.df_rna_tsne = read_tsne(self.match_dict["tsne_coord"])
        else:
            sys.stderr.write("rna tsne file not found!")

        # out
        self.mtx = f"{self.out_prefix}_citeseq.mtx.gz"
        self.matrix_dir = f"{self.out_prefix}_rna_citeseq_matrix"
        self.tsne_file = f"{self.out_prefix}_tsne_coord.tsv"

    @utils.add_log
    def run(self):
        mapped_read = int(self.df_read_count["read_count"].sum())

        # in cell
        df_read_count_in_cell = self.df_read_count[
            self.df_read_count.index.isin(self.match_barcode, level=0)
        ]
        mapped_read_in_cell = int(df_read_count_in_cell["read_count"].sum())
        self.add_metric(
            name="Mapped Reads in Cells",
            value=mapped_read_in_cell,
            total=mapped_read,
        )

        # UMI
        df_UMI_in_cell = df_read_count_in_cell.groupby(["barcode", TAG_COL]).agg(
            {"UMI": "count"}
        )

        df_temp = df_UMI_in_cell.reset_index().pivot(
            index="barcode", columns=TAG_COL, values="UMI"
        )
        df_cell = pd.DataFrame(index=self.match_barcode)
        df_pivot = pd.merge(
            df_cell, df_temp, how="left", left_index=True, right_index=True
        )

        # fillna
        df_pivot.fillna(0, inplace=True)
        df_pivot = df_pivot.astype(int)
        df_UMI_cell_out = df_pivot.T
        df_UMI_cell_out.to_csv(self.mtx, sep="\t", compression="gzip")

        # merge rna matrix
        tag_names = df_read_count_in_cell.index.get_level_values(1).unique()
        features = Features(tag_names)
        citeseq_matrix = CountMatrix.from_dataframe(
            df_UMI_in_cell, features, barcodes=self.match_barcode, value="UMI"
        )
        rna_matrix = CountMatrix.from_matrix_dir(matrix_dir=self.match_matrix_dir)
        merged_matrix = rna_matrix.concat_by_barcodes(citeseq_matrix)
        merged_matrix.to_matrix_dir(self.matrix_dir)

        # UMI
        UMIs = df_pivot.apply(sum, axis=1)
        median_umi = round(np.median(UMIs), 2)
        mean_umi = round(np.mean(UMIs), 2)
        self.add_metric(
            name="Median UMI per Cell",
            value=float(median_umi),
        )
        self.add_metric(
            name="Mean UMI per Cell",
            value=float(mean_umi),
        )

        # out tsne
        if not self.df_rna_tsne.empty:
            df_log1p = np.log2(df_pivot + 1)
            df_tsne = self.df_rna_tsne.merge(
                df_log1p, left_index=True, right_index=True
            )
            df_tsne.to_csv(self.tsne_file, sep="\t")


def get_opts_count_cite(parser, sub_program):
    if sub_program:
        parser.add_argument("--match_dir", help=HELP_DICT["match_dir"], required=True)
        parser.add_argument("--read_count_file", help="tag read count file")
        s_common(parser)


def count_cite(args):
    with Count_cite(args, display_title="Count") as runner:
        runner.run()
