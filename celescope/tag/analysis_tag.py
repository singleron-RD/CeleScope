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
    - Combine scRNA-Seq clustering infromation with tag assignment.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.tsne_tag_file = args.tsne_tag_file
        self.df_tsne_tag = pd.read_csv(self.tsne_tag_file, sep='\t')

    def run(self):
        report_runner = analysis_wrapper.Report_runner(self.args, display_title=self.display_title)
        df_tsne, _df_marker = report_runner.get_df()

        tsne_cluster = Tsne_plot(df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_tag = Tsne_plot(self.df_tsne_tag, 'tag').get_plotly_div()
        self.add_data(tsne_tag=tsne_tag)


def get_opts_analysis_tag(parser, sub_program):
    if sub_program:
        parser.add_argument('--tsne_tag_file', help='`{sample}_tsne_tag.tsv` from count_tag. ', required=True)
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--tsne_file", help=HELP_DICT['tsne_file'])
        parser = s_common(parser)


@utils.add_log
def analysis_tag(args):
    with Analysis_tag(args, display_title="Analysis") as runner:
        runner.run()
