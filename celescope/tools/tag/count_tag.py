from celescope.__init__ import HELP_DICT, HELP_INFO_DICT
from celescope.tools.step import Step, s_common
from celescope.tools import utils
import pandas as pd
import numpy as np


def get_opts_count_tag(parser, sub_program):
    if sub_program:
        parser.add_argument(
            "--read_count_file", help="Tag read count file.", required=True
        )
        parser.add_argument("--match_dir", help=HELP_DICT["match_dir"])
        parser.add_argument("--matrix_dir", help=HELP_DICT["matrix_dir"])
        parser.add_argument("--tsne_file", help=HELP_DICT["tsne_file"])

        s_common(parser)


class Count_tag(Step):
    """
    ## Features
    - Assign tag to each cell barcode and summarize.

    ## Output

    - `{sample}_umi_tag.tsv`

        `first column` cell barcode
        `last column`  assigned tag
        `columns between first and last` UMI count for each tag

    - `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.read_count_file = args.read_count_file

        # vals
        self.mapped_read = 0
        self.mapped_read_in_cell = 0
        self.df_UMI_cell = None
        self.df_tsne_tag = None

        # read
        self.df_read_count = pd.read_csv(self.read_count_file, sep="\t", index_col=0)

        if utils.check_arg_not_none(args, "match_dir"):
            match_dict = utils.parse_match_dir(args.match_dir)
            self.match_barcode = match_dict["match_barcode"]
            self.n_match_barcode = match_dict["n_match_barcode"]
            self.tsne_file = match_dict["tsne_coord"]
            self.matrix_dir = match_dict["matrix_dir"]
        elif utils.check_arg_not_none(args, "matrix_dir"):
            self.match_barcode, self.n_match_barcode = (
                utils.get_barcode_from_matrix_dir(args.matrix_dir)
            )
            self.tsne_file = args.tsne_file
            self.matrix_dir = args.matrix_dir
        else:
            raise ValueError("--match_dir or --matrix_dir is required.")

        # out files
        self.UMI_tag_file = f"{self.outdir}/{self.sample}_umi_tag.tsv"
        self.tsne_tag_file = f"{self.outdir}/{self.sample}_tsne_tag.tsv"

    def add_seq_metrics(self):
        self.add_metric(
            name=HELP_INFO_DICT["matched_barcode_number"]["display"],
            value=self.n_match_barcode,
            help_info=HELP_INFO_DICT["matched_barcode_number"]["info"],
        )

        self.add_metric(
            name="Mapped Reads in Cells",
            value=self.mapped_read_in_cell,
            total=self.mapped_read,
            help_info="Mapped reads with scRNA-Seq cell barcode",
        )

        UMIs = self.df_UMI_cell.apply(sum, axis=1)
        umi_median = round(np.median(UMIs), 2)
        umi_mean = round(np.mean(UMIs), 2)
        self.add_metric(
            name="Median UMI per Cell",
            value=umi_median,
            help_info="Median UMI per scRNA-Seq cell barcode",
        )
        self.add_metric(
            name="Mean UMI per Cell",
            value=umi_mean,
            help_info="Mean UMI per scRNA-Seq cell barcode",
        )

    def get_df_UMI_cell(self):
        self.mapped_read = int(self.df_read_count["read_count"].sum())
        df_read_count_in_cell = self.df_read_count[
            self.df_read_count.index.isin(self.match_barcode)
        ]
        self.mapped_read_in_cell = int(df_read_count_in_cell["read_count"].sum())
        tag_name = df_read_count_in_cell.columns[0]
        df_UMI_in_cell = (
            df_read_count_in_cell.reset_index()
            .groupby(["barcode", tag_name])
            .agg({"UMI": "count"})
        )
        df_UMI_in_cell = df_UMI_in_cell.reset_index()
        df_UMI_in_cell = df_UMI_in_cell.pivot(
            index="barcode", columns=tag_name, values="UMI"
        )
        df_cell = pd.DataFrame(index=self.match_barcode)
        df_UMI_cell = pd.merge(
            df_cell, df_UMI_in_cell, how="left", left_index=True, right_index=True
        )

        # fillna
        df_UMI_cell.fillna(0, inplace=True)
        self.df_UMI_cell = df_UMI_cell.astype(int)

    def get_df_tsne_tag(self):
        df_tsne = pd.read_csv(self.tsne_file, sep="\t", index_col=0)

        self.df_tsne_tag = pd.merge(
            df_tsne, self.df_UMI_cell, how="left", left_index=True, right_index=True
        )

    def write_files(self):
        self.df_UMI_cell.to_csv(self.UMI_tag_file, sep="\t")
        self.df_tsne_tag.to_csv(self.tsne_tag_file, sep="\t")

    def run(self):
        self.get_df_UMI_cell()
        self.get_df_tsne_tag()
        self.add_seq_metrics()
        self.write_files()


def count_tag(args):
    with Count_tag(args, display_title="Cells") as runner:
        runner.run()
