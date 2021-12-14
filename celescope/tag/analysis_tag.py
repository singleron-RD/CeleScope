import pandas as pd
import math
import plotly
import plotly.express as px
from collections import defaultdict
from celescope.snp.analysis_snp import Analysis_variant

import celescope.tools.utils as utils
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common


class Analysis_tag(AnalysisMixin):
    """
    Features
    - Combine scRNA-Seq clustering infromation with tag assignment.
    """

    def __init__(self, args,display_title=None):
        Step.__init__(self, args,display_title=display_title)
        AnalysisMixin.__init__(self, args)


    def plot_data(self,df):
        """
        Generate the data needed for the tsne_plot
        """
        tsne_data = df
        sum_cluster_df = tsne_data.groupby(['cluster']).agg("count").iloc[:, 0]
        sum_tag_df = tsne_data.groupby(['tag']).agg("count").iloc[:, 0]
        percent_cluster_df = sum_cluster_df.transform(lambda x: round(x / sum(x) * 100, 2))
        percent_tag_df = sum_tag_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_Clusters_dict = defaultdict(int)
        res_Clusters_list = []
        for cluster in sorted(tsne_data['cluster'].unique()):
            name = f"cluster {cluster}({percent_cluster_df[cluster]}%)"
            res_Clusters_dict[cluster]= name
            res_Clusters_list.append(name)
        tsne_data['name_Clusters'] = 'name'
        tsne_data['size'] = 4
        tsne_data.name_Clusters = tsne_data.cluster.map(res_Clusters_dict)
        res_tag_dict = defaultdict(int)
        res_tag_list = []
        for tag in sorted(tsne_data['tag'].unique()):
            name = f"{tag}({percent_tag_df[tag]}%)"
            res_tag_dict[tag]= name
            res_tag_list.append(name)    
        tsne_data['name_tag'] = 'name'
        tsne_data.name_tag = tsne_data.tag.map(res_tag_dict)
        return tsne_data,res_Clusters_list,res_tag_list

    @utils.add_log
    def type_plot(self,tsne_data,res_list,type):
        """
        Replace {type} of tsne_plot in HTML
        """
        x_ = math.ceil(max(abs(tsne_data["tSNE_1"])))
        y_ = math.ceil(max(abs(tsne_data["tSNE_2"])))
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
        fig = px.scatter(data_frame=tsne_data,
                                 title=f"t-SNE plot Colored by {type}",
                                 x="tSNE_1", 
                                 y="tSNE_2",
                                 color=f"name_{type}",
                                 size='size',
                                 size_max=4,
                                 labels={f"name_{type}":f"{type}"},
                                 category_orders={f"name_{type}":res_list})                               
        fig.update_yaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE2',
                         range=y_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig.update_xaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE1',
                         range=x_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig.update_layout(layout,title={"text":f"t-SNE plot Colored by {type}","x":0.5,"y":0.95,"font":{"size":15}},
                          plot_bgcolor = '#FFFFFF',hovermode="closest")
        div_item = plotly.offline.plot(fig, include_plotlyjs=True, output_type='div',config=config)
        return div_item

    @utils.add_log
    def add_plot(self,cluster,tag):
        downsample = {"downsample_clusters":cluster, "downsample_tag":tag}
        self.add_data(downsample=downsample)  


    def add_help(self):
        self.add_help_content(
            name='Marker Genes by Cluster',
            content='differential expression analysis based on the non-parameteric Wilcoxon rank sum test'
        )
        self.add_help_content(
            name='avg_log2FC',
            content='log fold-change of the average expression between the cluster and the rest of the sample'
        )
        self.add_help_content(
            name='pct.1',
            content='The percentage of cells where the gene is detected in the cluster'
        )
        self.add_help_content(
            name='pct.2',
            content='The percentage of cells where the gene is detected in the rest of the sample'
        )
        self.add_help_content(
            name='p_val_adj',
            content='Adjusted p-value, based on bonferroni correction using all genes in the dataset'
        )

    def run(self):
        self.get_analysis_data(feature_name="tag")
        tsne_tag_df = pd.read_csv(self.args.tsne_tag_file, sep="\t", index_col=0)
        feature_tsne = self.get_cluster_tsne(colname='tag', tsne_df=tsne_tag_df, show_colname=False)
        self.add_data(cluster_tsne=self.cluster_tsne)
        self.add_data(feature_tsne=feature_tsne)
        self.add_data(table_dict=self.table_dict)
        tsne_data,res_cluster_list,res_tag_list = self.plot_data(df=tsne_tag_df)
        div_cluster = self.type_plot(tsne_data,res_list=res_cluster_list,type="Clusters")
        div_tag = self.type_plot(tsne_data,res_list=res_tag_list,type="tag")
        self.add_plot(cluster=div_cluster,tag=div_tag)
        self.add_help()


def get_opts_analysis_tag(parser, sub_program):
    if sub_program:
        parser.add_argument('--tsne_tag_file', help='`{sample}_tsne_tag.tsv` from count_tag. ', required=True)
        parser.add_argument("--match_dir", help="Match celescope scRNA-Seq directory. ")
        parser.add_argument("--tsne_file", help="t-SNE coord file.")
        parser.add_argument("--marker_file", help="Marker file.")
        parser = s_common(parser)


@utils.add_log
def analysis_tag(args):
    with Analysis_tag(args,display_title="Analysis") as runner:
        runner.run()
