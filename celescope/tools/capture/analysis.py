import pandas as pd

from celescope.tools.step import s_common
from celescope.tools.plotly_plot import Tsne_plot
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.capture.__init__ import SUM_UMI_COLNAME



def get_opts_analysis(parser, sub_program):
    if sub_program:
        parser.add_argument('--filter_tsne_file', help='filter tsne file', required=True)
        s_common(parser)


class Analysis(AnalysisMixin):

    def __init__(self, args, display_title='Analysis'):
        super().__init__(args, display_title)
        self.df_tsne = pd.read_csv(args.filter_tsne_file)

    def add_cluster_metrics(self):
        self.add_help_content('cluster 1,2,3...', 'number of positive cells in each cluster after filtering')

        df_cluster_all = self.df_tsne.groupby("cluster").count()

        df_positive = self.df_tsne[self.df_tsne[SUM_UMI_COLNAME] > 0]
        df_cluster_positive = df_positive.groupby("cluster").count()

        for index, row in df_cluster_positive.iterrows():
            self.add_metric(
                name=f'cluster {index}',
                value=int(row['barcode']),
                total=int(df_cluster_all.loc[index, 'barcode']),
            )

    def run(self):
        self.add_cluster_metrics()

        tsne_cluster = Tsne_plot(self.df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_plot = Tsne_plot(self.df_tsne, SUM_UMI_COLNAME, discrete=False)
        tsne_plot.set_color_scale(['LightGrey', 'Orange', 'Red'])
        tsne_feature = tsne_plot.get_plotly_div()
        self.add_data(tsne_feature=tsne_feature)