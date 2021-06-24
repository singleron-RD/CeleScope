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

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        AnalysisMixin.__init__(self, args)
        self.matrix_file = args.matrix_file
        self.genomeDir = args.genomeDir
        self.type_marker_tsv = args.type_marker_tsv
        self.auto_assign_bool = False
        self.save_rds = args.save_rds
        if args.type_marker_tsv and args.type_marker_tsv != 'None':
            self.auto_assign_bool = True
            self.save_rds = True

    def run(self):
        self.seurat(self.matrix_file, self.save_rds, self.genomeDir)
        if self.auto_assign_bool:
            self.auto_assign(self.type_marker_tsv)
        self.run_analysis()
        self.add_data_item(cluster_tsne=self.cluster_tsne)
        self.add_data_item(gene_tsne=self.gene_tsne)
        self.add_data_item(table_dict=self.table_dict)

        self.clean_up()


@utils.add_log
def analysis(args):

    step_name = "analysis"
    ana = Analysis_rna(args, step_name)
    ana.run()


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
