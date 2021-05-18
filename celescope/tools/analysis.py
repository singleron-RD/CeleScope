import os
import subprocess

import pandas as pd

from celescope.tools.utils import add_log, gene_convert, s_common
from celescope.tools.analysis_mixin import AnalysisMixin 
from celescope.tools.step import Step


TOOLSDIR = os.path.dirname(__file__)


@add_log
def generate_matrix(gtf_file, matrix_file):

    id_name = gene_convert(gtf_file)
    matrix = pd.read_csv(matrix_file, sep="\t")

    gene_name_col = matrix.geneID.apply(lambda x: id_name[x])
    matrix.geneID = gene_name_col
    matrix = matrix.drop_duplicates(subset=["geneID"], keep="first")
    matrix = matrix.dropna()
    matrix = matrix.rename({"geneID": ""}, axis='columns')
    return matrix


@add_log
def seurat(sample, outdir, matrix_file, save_rds):
    app = TOOLSDIR + "/run_analysis.R"
    cmd = (
        f'Rscript {app} --sample {sample} --outdir {outdir} --matrix_file {matrix_file} '
        f'--save_rds {save_rds}'
    )
    seurat.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@add_log
def auto_assign(sample, outdir, type_marker_tsv):
    rds = f'{outdir}/{sample}.rds'
    app = TOOLSDIR + "/auto_assign.R"
    cmd = (
        f'Rscript {app} '
        f'--rds {rds} '
        f'--type_marker_tsv {type_marker_tsv} '
        f'--outdir {outdir} '
        f'--sample {sample} '
    )
    auto_assign.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


class Analysis_rna(Step, AnalysisMixin):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        AnalysisMixin.__init__(self, args)
        self.matrix_file = args.matrix_file
        self.type_marker_tsv = args.type_marker_tsv
        self.auto_assign_bool = False
        self.save_rds = args.save_rds
        if args.type_marker_tsv and args.type_marker_tsv != 'None':
            self.auto_assign_bool = True
            self.save_rds = True

    def run(self):
        seurat(self.sample, self.outdir, self.matrix_file, self.save_rds)
        if self.auto_assign_bool:
            auto_assign(self.sample, self.outdir, self.type_marker_tsv)
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
    parser.add_argument('--save_rds', action='store_true', help='write rds to disk')
    parser.add_argument('--type_marker_tsv', help='cell type marker tsv')

