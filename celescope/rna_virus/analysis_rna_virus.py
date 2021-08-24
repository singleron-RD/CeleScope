import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.step import Step
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.rna.analysis import get_opts_analysis


class Analysis_rna_virus(Step, AnalysisMixin):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        AnalysisMixin.__init__(self, args)

        # set
        self.virus_df = pd.read_csv(args.virus_file, sep="\t")

        self.auto_assign_bool = False
        if args.type_marker_tsv and args.type_marker_tsv != 'None':
            self.auto_assign_bool = True
            self.save_rds = True
    
    
    def run(self):
        self.seurat(self.args.matrix_file, self.args.save_rds, self.args.genomeDir)
        if self.auto_assign_bool:
            self.auto_assign(self.type_marker_tsv)
        self.run_analysis()
        virus_tsne = self.virus_tsne_list()
        self.add_data_item(cluster_tsne=self.cluster_tsne)
        self.add_data_item(virus_tsne=virus_tsne)
        self.add_data_item(table_dict=self.table_dict)

        self.clean_up()


    def virus_tsne_list(self):
        """
        return data dic
        """
        df = pd.merge(self.tsne_df, self.virus_df, on="barcode", how="left")
        df["UMI"] = df["UMI"].fillna(0)
        tSNE_1 = list(df.tSNE_1)
        tSNE_2 = list(df.tSNE_2)
        virus_UMI = list(df.UMI)
        res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "virus_UMI": virus_UMI}
        return res


@utils.add_log
def analysis_rna_virus(args):

    step_name = "analysis_rna_virus"
    runner = Analysis_rna_virus(args, step_name)
    runner.run()


def get_opts_analysis_rna_virus(parser, sub_program):
    get_opts_analysis(parser, sub_program)
    if sub_program:
        parser.add_argument(
            '--virus_file',
            help='virus UMI count file',
            required=True)
