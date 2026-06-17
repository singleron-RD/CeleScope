import numpy as np
import pandas as pd

from celescope.tools import utils, analysis_wrapper, plotly_plot
from celescope.tools.tag.analysis_tag import (
    get_opts_analysis_tag as opts,
)

# default colnames in tsne_tag file
DEFAULT_COLS = ["", "tSNE_1", "tSNE_2", "cluster", "Gene_Counts"]


class Analysis_tag(analysis_wrapper.Report_runner):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.umi_tag_file = args.umi_tag_file

        # out files
        self.tsne_tag_file = f"{self.outdir}/{self.sample}_tsne_tag.tsv"

    def run(self):
        df_tsne, _df_marker = self.get_df()

        # read umi_tag file and merge with tsne
        df_umi_tag = pd.read_csv(self.umi_tag_file, sep="\t", index_col=0)
        df_tsne_tag = pd.merge(
            df_tsne, df_umi_tag, how="left", left_index=True, right_index=True
        )

        # write tsne_tag file
        df_tsne_tag.to_csv(self.tsne_tag_file, sep="\t")

        colnames = list(df_tsne_tag.columns)
        tag_cols = set(colnames) - set(DEFAULT_COLS)
        if len(tag_cols) > 1:
            raise ValueError(
                f"Error: More than one tag barcode in tag_barcode_fasta: {tag_cols}. Currently, Multiple tag barcodes are not supported to display in HTML report."
            )
        tag_col = list(tag_cols)[0]

        log_tag_col = tag_col + "_log2"
        df_tsne_tag[log_tag_col] = np.log2(df_tsne_tag[tag_col] + 1)

        tsne_cluster = plotly_plot.Tsne_plot(df_tsne, "cluster").get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_tag = plotly_plot.Tsne_plot(
            df_tsne_tag, log_tag_col, discrete=False
        ).get_plotly_div()
        self.add_data(tsne_tag=tsne_tag)


@utils.add_log
def analysis_tag(args):
    with Analysis_tag(args, display_title="Analysis") as runner:
        runner.run()


def get_opts_analysis_tag(parser, sub_program):
    opts(parser, sub_program)
