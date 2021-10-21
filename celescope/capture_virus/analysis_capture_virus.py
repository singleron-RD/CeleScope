import os
import math

import numpy as np
import pandas as pd

import celescope.tools
from celescope.capture_virus.otsu import array2hist, makePlot, threshold_otsu
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common
import celescope.tools.utils as utils

TOOLS_DIR = os.path.dirname(celescope.tools.__file__)


@utils.add_log
def analysis_capture_virus(args):

    step_name = 'analysis_capture_virus'
    runner = Analysis_capture_virus(args, step_name)
    runner.run()


def get_opts_analysis_capture_virus(parser, sub_program):
    parser.add_argument("--umi_threshold", help='method to find virus UMI threshold',
                        choices=['otsu', 'auto', 'none'], default='otsu')
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='match_dir', required=True)
        parser.add_argument('--virus_file', help='virus UMI count file', required=True)


class Analysis_capture_virus(Step, AnalysisMixin):

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        AnalysisMixin.__init__(self, args)
        self.virus_file = args.virus_file
        self.umi_threshold = args.umi_threshold

        # parse
        _barcodes, self.n_cell = utils.read_barcode_file(self.match_dir)

        # set
        self.threshold = "None"
        self.add_metric(
            name='Number of cells',
            value=self.n_cell,
        )

        # out
        self.virus_tsne_file = f'{self.outdir}/{self.sample}_virus_tsne.tsv'
        self.auto_virus_file = f'{self.outdir}/{self.sample}_auto_UMI_count.tsv'
        self.otsu_virus_file = f'{self.outdir}/{self.sample}_otsu_UMI_count.tsv'
        self.otsu_plot = f'{self.outdir}/{self.sample}_otsu_plot.png'
    
    def run(self):
        if self.umi_threshold == 'otsu':
            self.otsu_threshold()
            virus_file = self.otsu_virus_file
        elif self.umi_threshold == 'auto':
            self.auto_threshold()
            virus_file = self.auto_virus_file
        else:
            virus_file = self.virus_file
        virus_df = pd.read_csv(virus_file, sep="\t")

        cluster_tsne = self.get_cluster_tsne(colname='cluster', tsne_df=self.tsne_df)
        virus_tsne = self.get_virus_tsne(virus_df)
        table_dict = self.get_table_dict(self.get_marker_table())

        self.add_data_item(cluster_tsne=cluster_tsne)
        self.add_data_item(virus_tsne=virus_tsne)
        self.add_data_item(table_dict=table_dict)
        self.clean_up()

    def get_virus_tsne(self, virus_df):
        virus_tsne_df = pd.merge(self.tsne_df, virus_df, on="barcode", how="left")
        virus_tsne_df.to_csv(self.virus_tsne_file, sep='\t')
        
        #get the sum of virus
        total_virus = virus_tsne_df.loc[:,"tag"].count()
        #get the sum of virus in each cluster
        cluster_virus = virus_tsne_df.groupby("cluster")["tag"].count()
        #Total cell each cluster
        total_cell_cluster = virus_tsne_df.groupby("cluster")["barcode"].count()
        total_cell = sum(total_cell_cluster)
        #add sum of virus to metric
        self.add_metric(name = 'Number of Cells with Virus',
                        value = total_virus,
                        total = total_cell)
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

        threshold = int(10 ** thresh)
        self.add_metric(
            name='Otsu UMI Threshold',
            value=threshold,
        )
        df_thresh = df_virus[df_virus["UMI"] >= threshold]
        df_thresh.to_csv(self.otsu_virus_file, sep='\t')

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



        

