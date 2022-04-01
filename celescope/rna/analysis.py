
from celescope.tools import analysis_wrapper
from celescope.tools.plotly_plot import Tsne_plot
from celescope.tools import utils
from celescope.tools.step import Step


class Analysis(Step):
    """
    ## Features
    - Cell clustering with Seurat.

    - Calculate the marker gene of each cluster.

    - Cell type annotation(optional). You can provide markers of known cell types and annotate cell types for each cluster.

    ## Output
    - `markers.tsv` Marker genes of each cluster.

    - `tsne_coord.tsv` t-SNE coordinates and clustering information.

    - `{sample}/06.analsis/{sample}_auto_assign/` This result will only be obtained when `--type_marker_tsv` 
    parameter is provided. The result contains 3 files:
        - `{sample}_auto_cluster_type.tsv` The cell type of each cluster; if cell_type is "NA", 
    it means that the given marker is not enough to identify the cluster.
        - `{sample}_png/{cluster}_pctdiff.png` Percentage of marker gene expression in this cluster - percentage in all other clusters.
        - `{sample}_png/{cluster}_logfc.png` log2 (average expression of marker gene in this cluster / average expression in all other clusters + 1)
    """

    def __init__(self, args, display_title=None):

        super().__init__(args, display_title)

        self.display_title = display_title
        self._table_id = 'marker_genes'

    def run(self):

        with analysis_wrapper.Scanpy_wrapper(self.args, display_title=self.display_title) as scanpy_wrapper:
            scanpy_wrapper.run()
            df_tsne, df_marker = scanpy_wrapper.get_df()
            self.set_metric_list(metric_list=scanpy_wrapper.get_metric_list())

        with analysis_wrapper.Report_runner(self.args, display_title=self.display_title) as report_runner:
            report_runner.add_marker_help()


        tsne_cluster = Tsne_plot(df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_gene = Tsne_plot(df_tsne, 'Gene_Counts', discrete=False).get_plotly_div()
        self.add_data(tsne_gene=tsne_gene)

        table_dict = self.get_table_dict(
            title='Marker Genes by Cluster',
            table_id=self._table_id,
            df_table=df_marker,
        )
        self.add_data(table_dict=table_dict)


@utils.add_log
def analysis(args):
    with Analysis(args, display_title='Analysis') as runner:
        runner.run()


def get_opts_analysis(parser, sub_program):
    analysis_wrapper.get_opts_analysis(parser, sub_program)
