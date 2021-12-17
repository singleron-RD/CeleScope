import pandas as pd
import plotly.express as px

from celescope.tools.step import Step, s_common
from celescope.tools.plotly_plot import Tsne_plot
from celescope.tools.analysis_mixin import AnalysisMixin
import celescope.tools.utils as utils


@utils.add_log
def analysis_virus(args):

    with Analysis_virus(args, display_title='Analysis') as runner:
        runner.run()

def get_opts_analysis_virus(parser, sub_program):
    if sub_program:
        parser.add_argument('--filter_tsne_file', help='filter tsne file', required=True)
        s_common(parser)


class Analysis_virus(AnalysisMixin):

    def __init__(self, args, display_title):
        super().__init__(args, display_title)
        self.df_tsne = pd.read_csv(args.filter_tsne_file)

    def run(self):
        tsne_cluster = Tsne_plot(self.df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_plot = Tsne_plot(self.df_tsne, 'UMI', discrete=False)
        tsne_plot.set_color_scale(['LightGrey', 'Orange', 'Red'])
        tsne_virus = tsne_plot.get_plotly_div()
        self.add_data(tsne_feature=tsne_virus)





        

