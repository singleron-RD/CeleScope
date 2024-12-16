import pandas as pd

from celescope.tools.capture.analysis import Analysis, get_opts_analysis
from celescope.fusion.count_fusion import Count_fusion
from celescope.fusion.mkref import Mkref_fusion
from celescope.tools.plotly_plot import Tsne_dropdown_plot, Tsne_plot


def analysis_fusion(args):
    with Analysis_fusion(args) as runner:
        runner.run()


def get_opts_analysis_fusion(parser, sub_program):
    parser.add_argument(
        "--fusion_genomeDir", help="Fusion genome directory.", required=True
    )
    get_opts_analysis(parser, sub_program)


class Analysis_fusion(Analysis):
    def __init__(self, args, display_title="Analysis"):
        super().__init__(args, display_title)

        fusion_pos_file = Mkref_fusion.get_config(args.fusion_genomeDir)["files"][
            "fusion_pos"
        ]
        self.pos_dict = Count_fusion.read_pos_file(fusion_pos_file)

        self.count_fusion_df = None
        self.count_fusion_out = f"{self.out_prefix}_fusion_count.csv"

    def plot_fusion(self):
        """
        plot fusion count
        """
        df_tsne_file = f"{self.out_prefix}_UMI_tsne.csv"
        df_tsne = pd.read_csv(df_tsne_file, sep=",", index_col=0)
        feature_name_list = df_tsne.columns[4:-1].to_list()
        # plot
        tsne_cluster = Tsne_plot(df_tsne, "cluster").get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)
        tsne_citeseq = Tsne_dropdown_plot(
            df_tsne, "Fusion", feature_name_list
        ).get_plotly_div()
        self.add_data(tsne_citeseq=tsne_citeseq)

    def get_fusion_count_df(self):
        fusion_df = self.df_tsne.reset_index().filter(
            list(self.pos_dict.keys()) + ["barcode"]
        )
        fusion_df = fusion_df.melt(
            id_vars=["barcode"], value_name="UMI", var_name="fusion", ignore_index=True
        )

        count = (
            fusion_df[fusion_df["UMI"] > 0]
            .groupby(["fusion"], as_index=False)["barcode"]
            .count()
        )
        self.count_fusion_df = count.assign(
            percent=lambda x: 100 * x.barcode / len(self.df_tsne)
        )
        self.count_fusion_df.to_csv(self.count_fusion_out, index=None)

    def add_fusion_count_metrics(self):
        self.add_help_content(
            "fusion_1, fusion_2...", "number of positive cells(UMI > 0) of each fusion"
        )
        for _, row in self.count_fusion_df.iterrows():
            self.add_metric(
                name=f'{row["fusion"]}',
                value=int(row["barcode"]),
                total=int(len(self.df_tsne)),
            )

    def run(self):
        super().run()
        self.plot_fusion()
        self.get_fusion_count_df()
        self.add_fusion_count_metrics()
