import numpy as np
import pandas as pd

from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.tools.plotly_plot import Tsne_dropdown_plot, Tsne_plot
from celescope.tools.utils import parse_match_dir, read_tsne


def get_opts_analysis_pathseq(parser, sub_program):
    if sub_program:
        parser.add_argument(
            "--umi_matrix_file", help="UMI matrix tsv file", required=True
        )
        parser.add_argument("--match_dir", help=HELP_DICT["match_dir"])
        parser.add_argument("--tsne_coord_file")
        s_common(parser)


def analysis_pathseq(args):
    with Analysis_pathseq(args, display_title="Analysis") as runner:
        runner.run()


class Analysis_pathseq(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        if args.match_dir:
            match_dict = parse_match_dir(self.args.match_dir)
            tsne_file = match_dict["tsne_coord"]
        else:
            tsne_file = args.tsne_coord_file
        self.df_tsne = read_tsne(tsne_file)
        self.df_umi = pd.read_csv(self.args.umi_matrix_file, sep="\t", index_col=0).T

    def add_dropdown_plot(self):
        sorted_columns = self.df_umi.sum().sort_values(ascending=False).index
        top10_feature = sorted_columns[:10]
        df_umi_log2 = self.df_umi.apply(lambda x: np.log2(x + 1), axis=0)
        df = self.df_tsne.merge(df_umi_log2, left_index=True, right_index=True)

        df["tSNE_1"] = df["tSNE_1"].astype("float16").round(2)
        df["tSNE_2"] = df["tSNE_2"].astype("float16").round(2)
        for col in top10_feature:
            df[col] = df[col].astype("float16").round(2)

        # plot
        tsne_cluster = Tsne_plot(df, "cluster").get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)
        tsne_feature = Tsne_dropdown_plot(
            df, "log2(UMI count+1)", top10_feature
        ).get_plotly_div()
        self.add_data(tsne_feature=tsne_feature)

    def add_summary_table(self):
        df = self.df_umi
        column_sums = df.sum()
        non_zero_counts = (df != 0).sum()
        # Create a new DataFrame with the information
        summary_df = pd.DataFrame(
            {
                "Genus": self.df_umi.columns,
                "Number of UMIs": column_sums,
                "Number of Positive Cells": non_zero_counts,
            }
        )
        summary_df = summary_df.sort_values(by="Number of UMIs", ascending=False)
        print(summary_df)
        table_dict = self.get_table_dict(
            title="Genus",
            table_id="pathseq",
            df_table=summary_df,
        )
        self.add_data(table_dict=table_dict)

    def run(self):
        self.add_summary_table()
        self.add_dropdown_plot()
