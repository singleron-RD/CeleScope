
from celescope.tools.analysis_mixin import AnalysisMixin, get_opts_analysis_mixin
from celescope.tools.plotly_plot import Tsne_plot
import celescope.tools.utils as utils


class Analysis(AnalysisMixin):
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

    def __init__(self, args, display_title=None):

        super().__init__(args, display_title)

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

        self._table_id = 'marker_genes'

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

        self.seurat(self.matrix_file, self.save_rds, self.genomeDir)

        self.add_mito_metric()
        if self.auto_assign_bool:
            self.auto_assign(self.type_marker_tsv)

        self.add_help()

        self.read_match_dir()

        tsne_cluster = Tsne_plot(self.df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_gene = Tsne_plot(self.df_tsne, 'Gene_Counts', discrete=False).get_plotly_div()
        self.add_data(tsne_gene=tsne_gene)

        table_dict = self.get_table_dict(
            title='Marker Genes by Cluster',
            table_id=self._table_id,
            df_table=self.df_marker,
        )
        self.add_data(table_dict=table_dict)


@utils.add_log
def analysis(args):
    with Analysis(args, display_title='Analysis') as runner:
        runner.run()


def get_opts_analysis(parser, sub_program):
    get_opts_analysis_mixin(parser, sub_program)
