from re import A
import pandas as pd

from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common
import celescope.tools.utils as utils


@utils.add_log
def generate_matrix(gtf_file, matrix_file):

    id_name = utils.get_id_name_dict(gtf_file)
    matrix = pd.read_csv(matrix_file, sep="\t")

    gene_name_col = matrix.geneID.apply(lambda x: id_name[x])
    matrix.geneID = gene_name_col
    matrix = matrix.drop_duplicates(subset=["geneID"], keep="first")
    matrix = matrix.dropna()
    matrix = matrix.rename({"geneID": ""}, axis='columns')
    return matrix


class Analysis_rna(Step, AnalysisMixin):
    """
    Features
    - Cell clustering with Seurat.

    - Calculate the marker gene of each cluster.

    - Cell type annotation(optional). You can provide markers of known cell types and annotate cell types for each cluster.

    Output
    - `markers.tsv` Marker genes of each cluster.

    - `tsne_coord.tsv` t-SNE coordinates and clustering information.

    - `{sample}/06.analsis/{sample}_auto_assign/` This result will only be obtained when `--type_marker_tsv` 
    parameter is provided. The result contains 3 files:
        - `{sample}_auto_cluster_type.tsv` The cell type of each cluster; if cell_type is "NA", 
    it means that the given marker is not enough to identify the cluster.
        - `{sample}_png/{cluster}_pctdiff.png` Percentage of marker gene expression in this cluster - percentage in all other clusters.
        - `{sample}_png/{cluster}_logfc.png` log2 (average expression of marker gene in this cluster / average expression in all other clusters + 1)
    """

    def __init__(self, args,display_title=None):
        Step.__init__(self, args, display_title=display_title)
        AnalysisMixin.__init__(self, args)
        self.matrix_file = args.matrix_file
        self.genomeDir = args.genomeDir
        self.type_marker_tsv = args.type_marker_tsv
        self.auto_assign_bool = False
        self.save_rds = args.save_rds
        self.outdir = args.outdir
        self.sample = args.sample
        if args.type_marker_tsv and args.type_marker_tsv != 'None':
            self.auto_assign_bool = True
            self.save_rds = True

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
        data = f'{self.outdir}/{self.sample}_tsne_coord.tsv'
        colnames = ["id","tsne_1","tsne_2","cluster","gene_counts"]
        tsne_data = pd.read_csv(data,sep="\t",header=0,names=colnames)
        sum_df = tsne_data.groupby(['cluster']).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_dict = defaultdict(int)
        res_list = []
        for cluster in sorted(tsne_data['cluster'].unique()):
            name = f"cluster {cluster}({percent_df[cluster]}%)"
            res_dict[cluster]= name
            res_list.append(name)
        tsne_data['name'] = 'name'
        tsne_data.name = tsne_data.cluster.map(res_dict)
        x_ = math.ceil(max(abs(tsne_data["tsne_1"])))
        y_ = math.ceil(max(abs(tsne_data["tsne_2"])))
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
        
        fig_clusters = px.scatter(data_frame=tsne_data,
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
        fig_clusters.update_layout(layout,title=dict(text= "t-SNE plot Colored by Clusters",x=0.4,y=0.95),
                          plot_bgcolor = '#FFFFFF',hovermode="closest")
        div_clusters = plotly.offline.plot(fig_clusters, include_plotlyjs=False, output_type='div',config=config)
        fig_gene = go.Figure()
        fig_gene.add_trace(go.Scatter(x = tsne_data['tsne_1'],y = tsne_data['tsne_2'],mode = 'markers',
                           marker = go.scatter.Marker(opacity=0.9,size=4,color=tsne_data['gene_counts'],colorscale='Jet',
                                                      colorbar=go.scatter.marker.ColorBar(title='Gene Counts')),
                           textposition='top center'))
        fig_gene.update_yaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE2',
                         range=y_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig_gene.update_xaxes(showgrid=True,gridcolor='#F5F5F5',showline=False, ticks=None,title_text='t-SNE1',
                         range=x_range,zeroline=True,zerolinecolor='black',zerolinewidth=0.7)
        fig_gene.update_layout(layout,title=dict(text= "t-SNE plot Colored by Gene Counts",x=0.5,y=0.95),
                          plot_bgcolor = '#FFFFFF',hovermode="closest")
        div_gene = plotly.offline.plot(fig_gene, include_plotlyjs=False, output_type='div',config=config)
        downsample = {"downsample_clusters":div_clusters, "downsample_gene":div_gene}
        self.add_data(downsample=downsample)        

    def run(self):

        self.seurat(self.matrix_file, self.save_rds, self.genomeDir)
        if self.auto_assign_bool:
            self.auto_assign(self.type_marker_tsv)

        self.get_analysis_data(feature_name="Gene Counts")
        self.add_data(cluster_tsne=self.cluster_tsne)
        self.add_data(feature_tsne=self.feature_tsne)
        self.add_data(table_dict=self.table_dict)
        self.add_help()
        self.tsne_plot()

@utils.add_log
def analysis(args):
    with Analysis_rna(args, display_title='Analysis') as runner:
        runner.run()


def get_opts_analysis(parser, sub_program):

    parser.add_argument('--genomeDir', help='Required. Genome directory.', required=True)
    parser.add_argument('--save_rds', action='store_true', help='Write rds to disk.')
    parser.add_argument(
        '--type_marker_tsv',
        help="""A tsv file with header. If this parameter is provided, cell type will be annotated. Example:
```
cell_type	marker
Alveolar	"CLDN18,FOLR1,AQP4,PEBP4"
Endothelial	"CLDN5,FLT1,CDH5,RAMP2"
Epithelial	"CAPS,TMEM190,PIFO,SNTN"
Fibroblast	"COL1A1,DCN,COL1A2,C1R"
B_cell	"CD79A,IGKC,IGLC3,IGHG3"
Myeloid	"LYZ,MARCO,FCGR3A"
T_cell	"CD3D,TRBC1,TRBC2,TRAC"
LUAD	"NKX2-1,NAPSA,EPCAM"
LUSC	"TP63,KRT5,KRT6A,KRT6B,EPCAM"
```"""
    )
    if sub_program:
        parser.add_argument(
            '--matrix_file',
            help='Required. Matrix_10X directory from step count.',
            required=True,
        )
        parser = s_common(parser)
