import pandas as pd

from celescope.tools.step import Step, s_common
from celescope.tools.plotly_plot import Tsne_dropdown_plot,Tsne_plot



def get_opts_analysis_cite(parser, sub_program):
    if sub_program:
        parser.add_argument('--tsne_coord', help='tsne coord file', required=True)
        s_common(parser)


def analysis_cite(args):
    with Analysis_cite(args, display_title='Analysis') as runner:
        runner.run()


class Analysis_cite(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)
        self.tmp_dir = self.args.outdir

        # input_file
        self.tsne_coord = self.args.tsne_coord

    def run(self):
        df_tsne = pd.read_csv(self.tsne_coord,sep="\t",index_col=0)
        feature_name_list = df_tsne.columns[4:].to_list()

        # accelerate
        df_tsne['tSNE_1'] = df_tsne['tSNE_1'].astype('float16')
        df_tsne['tSNE_2'] = df_tsne['tSNE_2'].astype('float16')
        for col in feature_name_list:
            df_tsne[col] = df_tsne[col].astype('float16')

        # plot
        tsne_cluster = Tsne_plot(df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)
        tsne_citeseq = Tsne_dropdown_plot(df_tsne,'Citeseq',feature_name_list,self.tmp_dir).get_plotly_div()
        self.add_data(tsne_citeseq=tsne_citeseq)
