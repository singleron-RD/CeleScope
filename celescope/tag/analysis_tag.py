import pandas as pd

from celescope.tools import utils
from celescope.tools.step import Step
from celescope.tools.step import s_common
from celescope.tools.plotly_plot import Tsne_plot
from celescope.tools import analysis_wrapper
from celescope.__init__ import HELP_DICT


class Analysis_tag(Step):
    """
    ## Features
    - Combine scRNA-Seq clustering information with tag assignment.

    ## Output

    - `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster information

    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.umi_tag_file = args.umi_tag_file

        # out files
        self.tsne_tag_file = f"{self.outdir}/{self.sample}_tsne_tag.tsv"

    def run(self):
        report_runner = analysis_wrapper.Report_runner(
            self.args, display_title=self.display_title
        )
        df_tsne, _df_marker = report_runner.get_df()

        tsne_cluster = Tsne_plot(df_tsne, "cluster").get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        # read umi_tag file and merge with tsne
        df_umi_tag = pd.read_csv(self.umi_tag_file, sep="\t", index_col=0)
        df_tsne_tag = pd.merge(
            df_tsne, df_umi_tag, how="left", left_index=True, right_index=True
        )

        # write tsne_tag file
        df_tsne_tag.to_csv(self.tsne_tag_file, sep="\t")

        tsne_tag = Tsne_plot(df_tsne_tag, "tag").get_plotly_div()
        self.add_data(tsne_tag=tsne_tag)


def get_opts_analysis_tag(parser, sub_program):
    if sub_program:
        parser.add_argument(
            "--umi_tag_file",
            help="`{sample}_umi_tag.tsv` from count_tag. ",
            required=True,
        )
        parser.add_argument("--match_dir", help=HELP_DICT["match_dir"])
        parser.add_argument("--tsne_file", help=HELP_DICT["tsne_file"])
        parser = s_common(parser)


@utils.add_log
def analysis_tag(args):
    with Analysis_tag(args, display_title="Analysis") as runner:
        runner.run()
