import pandas as pd

from celescope.tools.step import Step, s_common
from celescope.tools.plotly_plot import Tsne_dropdown_plot, Tsne_single_plot, Tsne_plot


def get_opts_analysis_pathseq(parser, sub_program):
    if sub_program:
        parser.add_argument(
            "--UMI_tsne", help="tsne coordinate plus UMI count file", required=True
        )
        s_common(parser)


def analysis_pathseq(args):
    with Analysis_pathseq(args, display_title="Analysis") as runner:
        runner.run()


class Analysis_pathseq(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)

    def run(self):
        df = pd.read_csv(self.args.UMI_tsne, sep="\t", index_col=0)
        feature_names = df.columns[4:]
        df_feature = df[feature_names]
        column_sums = df_feature.sum()
        non_zero_counts = (df_feature != 0).sum()
        # Create a new DataFrame with the information
        summary_df = pd.DataFrame(
            {"Number of UMIs": column_sums, "Number of Positive Cells": non_zero_counts}
        )
        summary_df = summary_df.sort_values(by="Number of UMIs", ascending=False)
        print(summary_df)
        sorted_columns = df_feature.sum().sort_values(ascending=False).index
        df_feature = df_feature[sorted_columns]

        # accelerate
        df["tSNE_1"] = df["tSNE_1"].astype("float16")
        df["tSNE_2"] = df["tSNE_2"].astype("float16")
        for col in feature_names:
            df[col] = df[col].astype("float16")

        # plot
        tsne_cluster = Tsne_plot(df, "cluster").get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)
        tsne_feature = Tsne_dropdown_plot(
            df, "UMI count", sorted_columns
        ).get_plotly_div()
        self.add_data(tsne_feature=tsne_feature)
        Tsne_single_plot(df, sorted_columns, self.args.outdir).get_plotly_div()
