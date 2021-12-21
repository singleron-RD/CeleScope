import pandas as pd

from celescope.tools.step import s_common
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

    def add_cluster_metrics(self):
        self.add_help_content('cluster 1,2,3...', 'number of cells with virus umi after filtering')

        df_cluster_all = self.df_tsne.groupby("cluster").count()

        df_virus = self.df_tsne[self.df_tsne['UMI'] > 0]
        df_cluster_virus = df_virus.groupby("cluster").count()

        for index, row in df_cluster_virus.iterrows():
            self.add_metric(
                name=f'cluster {index}',
                value=row['barcode'],
                total=df_cluster_all.loc[index, 'barcode'],
            )

    def run(self):
        self.add_cluster_metrics()

        tsne_cluster = Tsne_plot(self.df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_plot = Tsne_plot(self.df_tsne, 'UMI', discrete=False)
        tsne_plot.set_color_scale(['LightGrey', 'Orange', 'Red'])
        tsne_virus = tsne_plot.get_plotly_div()
        self.add_data(tsne_feature=tsne_virus)
