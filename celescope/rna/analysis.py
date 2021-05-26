import pandas as pd

from celescope.tools.utils import add_log, get_id_name_dict, s_common
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step



@add_log
def generate_matrix(gtf_file, matrix_file):

    id_name = get_id_name_dict(gtf_file)
    matrix = pd.read_csv(matrix_file, sep="\t")

    gene_name_col = matrix.geneID.apply(lambda x: id_name[x])
    matrix.geneID = gene_name_col
    matrix = matrix.drop_duplicates(subset=["geneID"], keep="first")
    matrix = matrix.dropna()
    matrix = matrix.rename({"geneID": ""}, axis='columns')
    return matrix


class Analysis_rna(Step, AnalysisMixin):
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


@add_log
def analysis(args):

    step_name = "analysis"
    ana = Analysis_rna(args, step_name)
    ana.run()


def get_opts_analysis(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--matrix_file', help='matrix file', required=True)
    parser.add_argument('--genomeDir', help='genomeDir', required=True)
    parser.add_argument('--save_rds', action='store_true', help='write rds to disk')
    parser.add_argument('--type_marker_tsv', help='cell type marker tsv')



