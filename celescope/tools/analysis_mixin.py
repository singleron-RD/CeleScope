import subprocess

import math
from collections import defaultdict

import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly.express as px

import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH
from celescope.rna.mkref import parse_genomeDir_rna
from celescope.tools.step import Step


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
            self.tsne_df = pd.read_csv(tsne_df_file, sep="\t")
            self.tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            marker_df_file = args.marker_file
            self.marker_df = pd.read_csv(marker_df_file, sep="\t")
        else:
            self.match_dir = args.outdir + "/../"  # use self


        # data
        self.table_dict = None
        self.feature_tsne = None

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

    @staticmethod
    def get_cluster_tsne(colname, tsne_df, show_colname=True):
        """
        tSNE_1	tSNE_2	cluster Gene_Counts
        return data list
        """

        sum_df = tsne_df.groupby([colname]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res = []
        for cluster in sorted(tsne_df[colname].unique()):
            sub_df = tsne_df[tsne_df[colname] == cluster]
            if show_colname:
                name = f"{colname} {cluster}({percent_df[cluster]}%)"
            else:
                name = f"{cluster}({percent_df[cluster]}%)"
            tSNE_1 = list(sub_df.tSNE_1)
            tSNE_2 = list(sub_df.tSNE_2)
            res.append({"name": name, "tSNE_1": tSNE_1, "tSNE_2": tSNE_2})
        return res

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

    def get_feature_tsne(self, feature_name):
        """
        return data dic
        """
        tsne_df = self.tsne_df
        tSNE_1 = list(tsne_df.tSNE_1)
        tSNE_2 = list(tsne_df.tSNE_2)
        Gene_Counts = list(tsne_df.Gene_Counts)
        res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "counts": Gene_Counts, 'feature_name':feature_name}
        self.feature_tsne = res

    def read_match_dir(self):
        """
        if match_dir is not self, should read match_dir at init
        if it is self, read at run_analysis - need to run seurat first
        """
        if self.match_dir:
            match_dict = utils.parse_match_dir(self.match_dir)
            tsne_df_file = match_dict['tsne_coord']
            self.tsne_df = pd.read_csv(tsne_df_file, sep="\t")
            self.tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            marker_df_file = match_dict['markers']
            self.marker_df = pd.read_csv(marker_df_file, sep="\t")

    @utils.add_log
    def continuous_tsne_plot(self):
        pass

    @utils.add_log
    def discrete_tsne_plot(self, feature_name):
        tsne_data = self.tsne_df
        sum_df = tsne_data.groupby([feature_name]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_dict = defaultdict(int)
        res_list = []
        for cluster in sorted(tsne_data[feature_name].unique()):
            name = f"{feature_name} {cluster}({percent_df[cluster]}%)"
            res_dict[cluster]= name
            res_list.append(name)
        tsne_data['name'] = 'name'
        tsne_data.name = tsne_data.cluster.map(res_dict)
        x_ = math.ceil(max(abs(tsne_data["tSNE_1"])))
        y_ = math.ceil(max(abs(tsne_data["tSNE_2"])))
        x_range,y_range = [-x_,x_],[-y_,y_]
        title=f"t-SNE plot Colored by {feature_name}"
        layout = {
            "height": 313,
            "width": 400,
            "margin": {
                "l": 45,
                "r": 35,
                "b": 30,
                "t": 30,}
                }
        config = {
            "displayModeBar": True, 
            "staticPlot": False, 
            "showAxisDragHandles": False, 
            "modeBarButtons": [["toImage", "resetScale2d"]], 
            "scrollZoom": True,
            "displaylogo": False
            }
        
        fig_clusters = px.scatter(data_frame=tsne_data,
                                 title=title,
                                 x="tSNE_1", 
                                 y="tSNE_2",
                                 color="name",
                                 labels={"name":feature_name},
                                 category_orders={"name":res_list})                               
        fig_clusters.update_yaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE2',
                         range=y_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig_clusters.update_xaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE1',
                         range=x_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig_clusters.update_layout(layout,title={"text":title, "x":0.5,"y":0.95,"font":{"size":15}},
                          plot_bgcolor = '#FFFFFF',hovermode="closest")
        div_clusters = plotly.offline.plot(fig_clusters, include_plotlyjs=True, output_type='div',config=config)

        return div_clusters


    @utils.add_log
    def tsne_plot(self):
        """
        Replace tsne_plot in HTML
        """
        import math
        import plotly
        import plotly.graph_objs as go
        import plotly.express as px
        from collections import defaultdict

        data = f'{self.out_prefix}_tsne_coord.tsv'
        colnames = ["id","tsne_1","tsne_2","cluster","gene_counts"]
        tsne_data = pd.read_csv(data,sep="\t",header=0,names=colnames)
        sum_df = tsne_data.groupby(['cluster']).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_dict = defaultdict(int)
        res_list = []
        for cluster in sorted(tsne_df['cluster'].unique()):
            name = f"cluster {cluster}({percent_df[cluster]}%)"
            res_dict[cluster]= name
            res_list.append(name)
        tsne_df['name'] = 'name'
        tsne_df.name = tsne_df.cluster.map(res_dict)
        x_ = math.ceil(max(abs(tsne_df["tsne_1"])))
        y_ = math.ceil(max(abs(tsne_df["tsne_2"])))
        x_range,y_range = [-x_,x_],[-y_,y_]
        layout = {
            "height": 313,
            "width": 400,
            "margin": {
                "l": 45,
                "r": 35,
                "b": 30,
                "t": 30,}
                }
        config = {
            "displayModeBar": True, 
            "staticPlot": False, 
            "showAxisDragHandles": False, 
            "modeBarButtons": [["toImage", "resetScale2d"]], 
            "scrollZoom": True,
            "displaylogo": False
            }
        
        fig_clusters = px.scatter(data_frame=tsne_df,
                                 title="t-SNE plot Colored by Clusters",
                                 x="tsne_1", 
                                 y="tsne_2",
                                 color="name",
                                 labels={"name":"cluster"},
                                 category_orders={"name":res_list})                               
        fig_clusters.update_yaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE2',
                         range=y_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig_clusters.update_xaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE1',
                         range=x_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig_clusters.update_layout(layout,title={"text":"t-SNE plot Colored by Clusters","x":0.5,"y":0.95,"font":{"size":15}},
                          plot_bgcolor = '#FFFFFF',hovermode="closest")
        div_clusters = plotly.offline.plot(fig_clusters, include_plotlyjs=True, output_type='div',config=config)
        fig_gene = go.Figure()
        fig_gene.add_trace(go.Scatter(x = tsne_df['tsne_1'],y = tsne_df['tsne_2'],mode = 'markers',
                           marker = go.scatter.Marker(opacity=0.9,size=4,color=tsne_df['gene_counts'],colorscale='Jet',
                                                      colorbar=go.scatter.marker.ColorBar(title='Gene Counts')),
                           textposition='top center'))
        fig_gene.update_yaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE2',
                         range=y_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig_gene.update_xaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE1',
                         range=x_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig_gene.update_layout(layout,title={"text":"t-SNE plot Colored by Gene Counts","x":0.5,"y":0.95,"font":{"size":15}},
                          plot_bgcolor = '#FFFFFF',hovermode="closest")
        div_gene = plotly.offline.plot(fig_gene, include_plotlyjs=True, output_type='div',config=config)
        downsample = {"downsample_clusters":div_clusters, "downsample_gene":div_gene}
        self.add_data(downsample=downsample)   

    def get_analysis_data(self, feature_name):
        self.read_match_dir()
        self.cluster_tsne = self.get_cluster_tsne(colname='cluster', tsne_df=self.tsne_df)
        self.get_feature_tsne(feature_name)
        marker_table = self.get_marker_table()
        self.get_table_dict(marker_table)
