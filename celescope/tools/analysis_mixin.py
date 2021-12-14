import subprocess

import math
from collections import defaultdict

import pandas as pd
import plotly
import plotly.express as px

import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH
from celescope.rna.mkref import parse_genomeDir_rna
from celescope.tools.step import Step


PLOTLY_CONFIG =  {
    "displayModeBar": True, 
    "staticPlot": False, 
    "showAxisDragHandles": False, 
    "modeBarButtons": [["toImage", "resetScale2d"]], 
    "scrollZoom": False,
    "displaylogo": False
}


class Tsne_plot():

    def __init__(self, df_tsne, feature_name, discrete=True):
        self.df_tsne = df_tsne
        self.feature_name = feature_name
        self.discrete = discrete
        self.title = f"t-SNE plot Colored by {feature_name}"
        
        self._layout = {}
        self._dot_size = 4
        self._str_coord1 = "tSNE_1"
        self._str_coord2 = "tSNE_2"
        self.axes_config = {
            'showgrid': True,
            'gridcolor': '#F5F5F5',
            'showline': False, 
            'ticks': None,
            'zeroline': True,
            'zerolinecolor': 'black',
            'zerolinewidth': 0.7,
        }
        x_ = math.ceil(max(abs(df_tsne[self._str_coord1])))
        y_ = math.ceil(max(abs(df_tsne[self._str_coord2])))
        self.x_range, self.y_range = [-x_,x_], [-y_,y_] 

        self.scatter_config = {
            'data_frame': df_tsne,
            'title': self.title,
            'x': self._str_coord1, 
            'y': self._str_coord2,
            'size_max': self._dot_size,
            'opacity': 0.9,
            'color': self.feature_name,
            'color_continuous_scale': px.colors.sequential.Jet,
        }

        self._fig = None

        if discrete:
            self.discrete_tsne_plot()
        else:
            self.continuous_tsne_plot()
        self.update_fig()

        self.plotly_div = plotly.offline.plot(self._fig, include_plotlyjs=True, output_type='div', config=PLOTLY_CONFIG)

    @utils.add_log
    def discrete_tsne_plot(self):
        
        df_tsne = self.df_tsne
        feature_name = self.feature_name

        sum_df = df_tsne.groupby([feature_name]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_dict = defaultdict(int)
        res_list = []
        for cluster in sorted(df_tsne[feature_name].unique()):
            name = f"{cluster}({percent_df[cluster]}%)"
            res_dict[cluster]= name
            res_list.append(name)

        df_tsne[self.feature_name] = df_tsne[self.feature_name].map(res_dict)

        self._fig = px.scatter(
            **self.scatter_config,
            category_orders={"name":res_list}
        )                               

    @utils.add_log
    def continuous_tsne_plot(self):

        self._fig = px.scatter(
            **self.scatter_config,
        )   



    def update_fig(self):
        self._fig.update_xaxes(
            title_text=self._str_coord1,
            range=self.x_range,
            **self.axes_config
        )

        self._fig.update_yaxes(
            title_text=self._str_coord2,
            range=self.y_range,
            **self.axes_config
        )
        
        self._fig.update_layout(
            self._layout,
            title={ "text":self.title, "x":0.5, "y":0.95, "font":{"size":15} },
            plot_bgcolor = '#FFFFFF',
            hovermode="closest"
        )

class AnalysisMixin(Step):
    """
    mixin class for analysis
    """

    def __init__(self, args, display_title=None):
        
        super().__init__(args, display_title=display_title)

        self.match_dir = None
        if hasattr(args, "match_dir") and args.match_dir:
            self.match_dir = args.match_dir
            self.read_match_dir()
        elif hasattr(args, "tsne_file") and args.tsne_file:
            tsne_df_file = args.tsne_file
            self.df_tsne = pd.read_csv(tsne_df_file, sep="\t")
            self.df_tsne.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            marker_df_file = args.marker_file
            self.marker_df = pd.read_csv(marker_df_file, sep="\t")
        else:
            self.match_dir = args.outdir + "/../"  # use self


        # data
        self.table_dict = None

    @utils.add_log
    def seurat(self, matrix_file, save_rds, genomeDir):
        app = ROOT_PATH + "/tools/run_analysis.R"
        genome = parse_genomeDir_rna(genomeDir)
        mt_gene_list = genome['mt_gene_list']
        cmd = (
            f'Rscript {app} '
            f'--sample {self.sample} '
            f'--outdir {self.outdir} '
            f'--matrix_file {matrix_file} '
            f'--mt_gene_list {mt_gene_list} '
            f'--save_rds {save_rds}'
        )
        self.seurat.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def auto_assign(self, type_marker_tsv):
        rds = f'{self.outdir}/{self.sample}.rds'
        app = ROOT_PATH + "/tools/auto_assign.R"
        cmd = (
            f'Rscript {app} '
            f'--rds {rds} '
            f'--type_marker_tsv {type_marker_tsv} '
            f'--outdir {self.outdir} '
            f'--sample {self.sample} '
        )
        self.auto_assign.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)


    def get_table_dict(self, marker_table):
        self.table_dict = Step.get_table(
            title='Marker Genes by Cluster',
            table_id='marker_gene_table',
            df_table=marker_table,
        )

    def get_marker_table(self):
        """
        return html code
        """

        avg_logfc_col = "avg_log2FC"  # seurat 4
        if "avg_logFC" in self.marker_df.columns:  # seurat 2.3.4
            avg_logfc_col = "avg_logFC"
        marker_df = self.marker_df.loc[:,
                                       ["cluster", "gene", avg_logfc_col, "pct.1", "pct.2", "p_val_adj"]
                                       ]
        marker_df["cluster"] = marker_df["cluster"].apply(lambda x: f"cluster {x}")

        return marker_df


    def read_match_dir(self):
        """
        if match_dir is not self, should read match_dir at init
        if it is self, read at run_analysis - need to run seurat first
        """
        if self.match_dir:
            match_dict = utils.parse_match_dir(self.match_dir)
            tsne_df_file = match_dict['tsne_coord']
            self.df_tsne = pd.read_csv(tsne_df_file, sep="\t")
            self.df_tsne.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            marker_df_file = match_dict['markers']
            self.marker_df = pd.read_csv(marker_df_file, sep="\t")


