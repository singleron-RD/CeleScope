import os
import math

import numpy as np
import pandas as pd
from celescope.rna.analysis import Analysis

import celescope.tools
from celescope.capture_virus.otsu import array2hist, makePlot, threshold_otsu
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common
import celescope.tools.utils as utils

TOOLS_DIR = os.path.dirname(celescope.tools.__file__)


@utils.add_log
def analysis_capture_virus(args):

    with Analysis_capture_virus(args) as runner:
        runner.run()

def get_opts_analysis_capture_virus(parser, sub_program):
    parser.add_argument("--umi_threshold_method", help='method to find virus UMI threshold',
                        choices=['otsu', 'auto', 'none'], default='otsu')
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='match_dir', required=True)
        parser.add_argument('--virus_file', help='virus UMI count file', required=True)


class Analysis_capture_virus(AnalysisMixin):

    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        self.virus_file = args.virus_file
        self.df_virus= pd.read_csv(virus_file, sep="\t")
        self.umi_threshold_method = args.umi_threshold_method

        # parse
        _barcodes, self.n_cell = utils.read_barcode_file(self.match_dir)

        # set
        self.umi_threshold = 1

        # out
        self.filter_virus_tsne_file = f'{self.outdir}/{self.sample}_filtered_virus_tsne.tsv'
        self.otsu_plot = f'{self.outdir}/{self.sample}_otsu_plot.png'
    
    def get_umi_threshold(self):
        if self.umi_threshold_method == 'auto':
            self.auto_threshold()
        elif self.umi_threshold_method == 'otsu':
            self.otsu_threshold()


    def filter_write_filter_tsne(self):
        df_filter = self.df_virus[self.df_virus["UMI"] >= threshold]
        df_filter_tsne = pd.merge(self.df_tsne, virus_df, on="barcode", how="left")
        df_filter.fillna(0, inplace=True)
        df_filter.to_csv(self.auto_virus_file, sep='\t')

    def run(self):
        self.get_umi_threshold()
        self.filter_write_virus_df()

        self.add_data(tsne_cluster=tsne_cluster)
        self.add_data(tsne_virus=tsne_virus)


    def get_virus_tsne(self, virus_df):
        virus_tsne_df = pd.merge(self.df_tsne, virus_df, on="barcode", how="left")
        virus_tsne_df.to_csv(self.virus_tsne_file, sep='\t')
        
        #get the sum of virus
        total_virus = virus_tsne_df.loc[:,"tag"].count()
        #get the sum of virus in each cluster
        cluster_virus = virus_tsne_df.groupby("cluster")["tag"].count()
        #Total cell each cluster
        total_cell_cluster = virus_tsne_df.groupby("cluster")["barcode"].count()
        total_cell = sum(total_cell_cluster)

        self.add_metric(
            name='Number of Cells with Virus',
            value=total_virus,
            total=total_cell,
            help_info='Number of cells detected with virus',
        )
        cluster = 1    
        for num_virus in cluster_virus:
            cluster_cell_num = total_cell_cluster[cluster]
            self.add_metric(name = f'Cluster_{cluster}',
                            value = num_virus,
                            total = cluster_cell_num)
            cluster += 1


        virus_tsne_df["UMI"] = virus_tsne_df["UMI"].fillna(0)
        tSNE_1 = list(virus_tsne_df.tSNE_1)
        tSNE_2 = list(virus_tsne_df.tSNE_2)
        virus_UMI = list(virus_tsne_df.UMI)
        res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "virus_UMI": virus_UMI}
        return res

    @utils.add_log
    def otsu_threshold(self):
        df_virus = pd.read_csv(self.virus_file, sep='\t')
        array = np.log10(df_virus["UMI"])
        hist = array2hist(array)
        thresh = threshold_otsu(hist)
        makePlot(hist, thresh, self.otsu_plot)

        self.threshold = int(10 ** thresh)

    @utils.add_log
    def auto_threshold(self):
        """
        threhold = 99 percentile of all cell virus UMIs / 10
        """
        umis = self.df_virus["UMI"]

        cell_99th = len(umis) // 100
        sorted_umis = sorted(umis, reverse=True)
        percentile_99_umi = sorted_umis[cell_99th]
        self.threshold = int(percentile_99_umi / 10)
    

    def add_metrics(self):

        self.add_metric(
            name='Number of cells',
            value=self.n_cell,
            help_info='Number of cells in the matched scRNA-Seq sample',
        )

        self.add_metric(
            name='UMI Threshold Method',
            value=self.umi_threshold_method,
            help_info='Method to find UMI threshold',
        )

        self.add_metric(
            name='UMI Threshold',
            value=self.umi_threshold,
            help_info='UMI threshold',
        )
        df_threshold = df_virus[df_virus["UMI"] >= threshold]
        df_threshold.to_csv(self.otsu_virus_file, sep='\t')



    @utils.add_log
    def auto_threshold(self):
        """
        threhold = 99 percentile of all cell virus UMIs / 10
        """
        df_virus = pd.read_csv(self.virus_file, sep='\t')
        umis = df_virus["UMI"]

        cell_99th = self.n_cell // 100
        if cell_99th >= len(umis):
            threshold = 1
        else:
            sorted_umis = sorted(umis, reverse=True)
            percentile_99_umi = sorted_umis[cell_99th]
            threshold = math.ceil(percentile_99_umi / 10)
        self.add_metric(
            name='Auto UMI Threshold',
            value=threshold,
        )
        df_threshold = df_virus[df_virus["UMI"] >= threshold]
        df_threshold.to_csv(self.auto_virus_file, sep='\t')



        

