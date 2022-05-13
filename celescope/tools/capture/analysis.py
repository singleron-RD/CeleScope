import pandas as pd

from celescope.tools.plotly_plot import Tsne_plot
from celescope.tools.step import Step
from celescope.tools.capture.__init__ import SUM_UMI_COLNAME
from celescope.tools import analysis_wrapper



def get_opts_analysis(parser, sub_program):
    if sub_program:
        parser.add_argument('--filter_umi_file', help='filter umi file', required=True)
        analysis_wrapper.get_opts_analysis_match(parser, sub_program)


class Analysis(Step):

    def __init__(self, args, display_title='Analysis'):
        super().__init__(args, display_title)
        self.df_filter_umi = pd.read_csv(args.filter_umi_file, index_col=0)

        # out
        self.df_tsne_file = f'{self.out_prefix}_UMI_tsne.csv'
        self.df_tsne = None

    def add_tsne_info(self):

        report_runner = analysis_wrapper.Report_runner(self.args)
        df_tsne, _df_marker = report_runner.get_df()
        self.df_tsne = pd.merge(df_tsne, self.df_filter_umi, left_index=True, right_index=True)
        self.df_tsne.to_csv(self.df_tsne_file)

    def add_cluster_metrics(self):
        self.add_help_content('cluster 1,2,3...', 'number of positive cells(sum_UMI > 0) in each cluster after filtering')

        df_cluster_all = self.df_tsne.groupby("cluster").count()

        df_positive = self.df_tsne[self.df_tsne[SUM_UMI_COLNAME] > 0]
        df_cluster_positive = df_positive.groupby("cluster").count()

        for index, row in df_cluster_positive.iterrows():
            self.add_metric(
                name=f'cluster {index}',
                value=int(row['tSNE_1']),
                total=int(df_cluster_all.loc[index, 'tSNE_1']),
            )

    def run(self):
        self.add_tsne_info()
        self.add_cluster_metrics()

        tsne_cluster = Tsne_plot(self.df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_plot = Tsne_plot(self.df_tsne, SUM_UMI_COLNAME, discrete=False)
        tsne_plot.set_color_scale(['LightGrey', 'Orange', 'Red'])
        tsne_feature = tsne_plot.get_plotly_div()
        self.add_data(tsne_feature=tsne_feature)