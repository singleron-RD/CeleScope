import pandas as pd

from celescope.tools import utils
from celescope.tools.step import Step
from celescope.rna.analysis import get_opts_analysis


class Analysis_rna_virus(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)

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
        self.get_analysis_data()
        virus_tsne = self.virus_tsne_list()
        self.add_data(cluster_tsne=self.cluster_tsne)
        self.add_data(virus_tsne=virus_tsne)
        self.add_data(table_dict=self.table_dict)

        self._clean_up()

    def virus_tsne_list(self):
        """
        return data dic
        """
        df = pd.merge(self.df_tsne, self.virus_df, on="barcode", how="left")
        df["UMI"] = df["UMI"].fillna(0)
        tSNE_1 = list(df.tSNE_1)
        tSNE_2 = list(df.tSNE_2)
        virus_UMI = list(df.UMI)
        res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "virus_UMI": virus_UMI}
        return res


@utils.add_log
def analysis_rna_virus(args):

    with Analysis_rna_virus(args) as runner:
        runner.run()


def get_opts_analysis_rna_virus(parser, sub_program):
    get_opts_analysis(parser, sub_program)
    if sub_program:
        parser.add_argument(
            '--virus_file',
            help='virus UMI count file',
            required=True)
