import numpy as np
import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common
import celescope.capture_virus.otsu as otsu 
from celescope.tools.plotly_plot import Tsne_plot


def get_opts_filter_virus(parser, sub_program):

    parser.add_argument(
        "--umi_threshold_method", 
        help='method to find virus UMI threshold',
        choices=['otsu', 'auto', 'hard'], 
        default='otsu'
    )
    parser.add_argument(
        "--umi_hard_threshold", 
        help='int, use together with `--umi_threshold_method hard`',
    )
    if sub_program:
        parser.add_argument('--umi_tsne_file', required=True)
        s_common(parser)

def filter_virus(args):
    
    with Filter_virus(args, display_title='Filtering') as runner:
        runner.run()


class Filter_virus(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)
        self.umi_threshold_method = args.umi_threshold_method
        self.df_umi_tsne = pd.read_csv(args.umi_tsne_file)
        
        # data
        self.threshold = 1 # if not set explicitly, use 1 as default
        self.df_filter = self.df_umi_tsne.fillna({'UMI': 0})
        self.umi_array = [x for x in self.df_filter['UMI'].values if x > 0]

        # out
        self.otsu_plot = f'{self.out_prefix}_otsu.png'
        self.filter_tsne_file = f'{self.out_prefix}_filter_tsne.csv'

    def get_umi_threshold(self):
        if self.umi_threshold_method == 'auto':
            self.auto_threshold()
        elif self.umi_threshold_method == 'otsu':
            self.otsu_threshold()
        elif self.umi_threshold_method == 'hard':
            self.hard_threshold()

        self.add_metric('UMI threshold method', self.umi_threshold_method)
        self.add_metric('UMI threshold', self.threshold)

    @utils.add_log
    def otsu_threshold(self):
        array = np.log2(self.umi_array)
        hist = otsu.array2hist(array)
        thresh = otsu.threshold_otsu(hist)
        otsu.makePlot(hist, thresh, self.otsu_plot)

        self.threshold = int(2 ** thresh)

    @utils.add_log
    def auto_threshold(self):
        """
        threhold = 99 percentile of all cell virus UMIs / 10
        """
        umi_array = self.umi_array

        cell_99th = len(umi_array) // 100
        sorted_umis = sorted(umi_array, reverse=True)
        percentile_99_umi = sorted_umis[cell_99th]
        self.threshold = int(percentile_99_umi / 10)
    
    @utils.add_log
    def hard_threshold(self):
        self.threshold = self.args.umi_hard_threshold


    def filter_write(self):
        self.df_filter.loc[self.df_filter["UMI"] < self.threshold, 'UMI'] = 0
        self.df_filter.to_csv(self.filter_tsne_file)

        cell_total = len(self.df_filter)
        df_virus = self.df_filter[self.df_filter['UMI'] > 0]
        cell_with_virus = len(df_virus)
        self.add_metric(
            name='Number of Cells with Virus',
            value=cell_with_virus,
            total=cell_total,
        )

        df_cluster_all = self.df_filter.groupby("cluster").count()
        df_cluster_virus = df_virus.groupby("cluster").count()
        for index, row in df_cluster_virus.iterrows():
            self.add_metric(
                name=f'cluster {index}',
                value=row['barcode'],
                total=df_cluster_all.loc[index, 'barcode'],
            )



    def run(self):
        self.get_umi_threshold()
        self.filter_write()


