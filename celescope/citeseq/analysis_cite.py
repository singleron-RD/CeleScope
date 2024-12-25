import pandas as pd

from celescope.tools.step import Step, s_common
from celescope.tools.plotly_plot import Tsne_dropdown_plot, Tsne_plot


def get_opts_analysis_cite(parser, sub_program):
    if sub_program:
        parser.add_argument("--tsne_coord", help="tsne coord file", required=True)
        s_common(parser)


def analysis_cite(args):
    with Analysis_cite(args, display_title="Analysis") as runner:
        runner.run()


class Analysis_cite(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)

    def run(self):
        df_tsne = pd.read_csv(self.args.tsne_coord, sep="\t", index_col=0)
        # filter low expression
        keep_antibody = df_tsne.columns[4:][
            (df_tsne.iloc[:, 4:].mean() > 0.1).to_list()
        ].to_list()
        keep_cols = df_tsne.columns[:4].to_list() + keep_antibody
        df_tsne = df_tsne.loc[:, keep_cols]

        # accelerate
        df_tsne["tSNE_1"] = df_tsne["tSNE_1"].astype("float16")
        df_tsne["tSNE_2"] = df_tsne["tSNE_2"].astype("float16")
        for col in keep_antibody:
            df_tsne[col] = df_tsne[col].astype("float16")

        # plot
        tsne_cluster = Tsne_plot(df_tsne, "cluster").get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)
        tsne_citeseq = Tsne_dropdown_plot(
            df_tsne, "Citeseq", keep_antibody
        ).get_plotly_div()
        self.add_data(tsne_citeseq=tsne_citeseq)
