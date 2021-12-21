import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import s_common
from celescope.tools.plotly_plot import Tsne_plot


class Analysis_tag(AnalysisMixin):
    """
    Features
    - Combine scRNA-Seq clustering infromation with tag assignment.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.tsne_tag_file = args.tsne_tag_file
        self.df_tsne_tag = pd.read_csv(self.tsne_tag_file, sep='\t')

    def run(self):

        tsne_cluster = Tsne_plot(self.df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_tag = Tsne_plot(self.df_tsne_tag, 'tag').get_plotly_div()
        self.add_data(tsne_tag=tsne_tag)


def get_opts_analysis_tag(parser, sub_program):
    if sub_program:
        parser.add_argument('--tsne_tag_file', help='`{sample}_tsne_tag.tsv` from count_tag. ', required=True)
        parser.add_argument("--match_dir", help="Match celescope scRNA-Seq directory. ")
        parser.add_argument("--tsne_file", help="t-SNE coord file.")
        parser = s_common(parser)


@utils.add_log
def analysis_tag(args):
    with Analysis_tag(args, display_title="Analysis") as runner:
        runner.run()
