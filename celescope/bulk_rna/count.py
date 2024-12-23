import random
import pandas as pd
import numpy as np
from celescope.tools import utils
from celescope.tools import count as super_count
from celescope.tools.matrix import ROW
from celescope.tools.plotly_plot import Line_plot

random.seed(0)
np.random.seed(0)


class Count(super_count.Count):
    """
    ## Features
    - Generate expression matrix.
    - Well statstic.

    ## Output
    - `{sample}_raw_feature_bc_matrix` The expression matrix of all detected barcodes in [Matrix Market Exchange Formats](
        https://math.nist.gov/MatrixMarket/formats.html).
    - `{sample}_count_detail.txt.gz` 4 columns:
        - barcode
        - gene ID
        - UMI count
        - read_count
    - `{sample}_counts.txt` 6 columns:
        - Barcode: barcode sequence
        - read: read count of each barcode
        - UMI: UMI count for each barcode
        - geneID: gene count for each barcode
        - mark: cell barcode or backgound barcode.
            `CB` cell
            `UB` background
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        self.raw_count_file = f"{self.outdir}/{self.sample}_counts.txt"
        self.marked_count_file = f"{self.outdir}/{self.sample}_counts_report.txt"
        self.umi_cutoff = args.umi_cutoff
        self.read_cutoff = args.read_cutoff
        self.gene_cutoff = args.gene_cutoff
        self._table_id = self.assay

    @utils.add_log
    def run(self):
        ## output exprssion matrix
        df = pd.read_table(
            self.args.count_detail, header=0, dtype={"geneID": str}, index_col=[0, 1]
        )
        self.write_sparse_matrix(df, self.raw_matrix_dir)
        ## output stats
        df_bc = self.get_df_bc(df)
        df_bc.index.name = "Well"
        sort_col = ["UMI", "read", ROW]
        df_bc = df_bc[sort_col]
        df_bc.columns = ["UMI", "read", "gene"]
        df_bc.to_csv(self.raw_count_file, sep="\t")
        df_bc_valid = df_bc[
            (df_bc["UMI"] >= self.umi_cutoff)
            & (df_bc["read"] >= self.read_cutoff)
            & (df_bc["gene"] >= self.gene_cutoff)
        ]
        df_bc_valid.to_csv(self.marked_count_file, sep="\t")

        ## saturation
        df_valid = (
            df.reset_index().set_index("Barcode").loc[df_bc_valid.index.to_list()]
        )
        self.saturation = Count.get_read_saturation(df_valid)
        self.downsample_dict = self.downsample(df_valid.reset_index())

        ## report
        self.add_count_metrics(df_bc_valid)
        self.add_plot_data()

    @staticmethod
    def get_read_saturation(df_cell):
        unique = sum(df_cell["unique"])
        reads = sum(df_cell["read"])
        saturation = 1 - unique / reads
        return saturation

    @staticmethod
    def sub_sample(fraction, df_cell, cell_read_index):
        cell_read = df_cell["read"].sum()
        frac_n_read = int(cell_read * fraction)
        subsample_read_index = cell_read_index[:frac_n_read]
        index_dedup = np.unique(subsample_read_index, return_counts=False)
        # gene median
        df_cell_subsample = df_cell.loc[index_dedup,]
        geneNum_median = float(
            df_cell_subsample.groupby("Barcode").agg({ROW: "nunique"}).median().iloc[0]
        )
        return geneNum_median

    @utils.add_log
    def downsample(self, df_cell):
        """saturation and median gene
        return fraction=1 saturation
        """
        cell_read_index = np.array(df_cell.index.repeat(df_cell["read"]), dtype="int32")
        np.random.shuffle(cell_read_index)

        format_str = "%.2f\t%.2f\n"
        res_dict = {"fraction": [0.0], "median_gene": [0]}
        with open(self.downsample_file, "w") as fh:
            fh.write("percent\tmedian_geneNum\n")
            fh.write(format_str % (0, 0))
            for fraction in np.arange(0.1, 1.1, 0.1):
                geneNum_median = Count.sub_sample(fraction, df_cell, cell_read_index)
                fh.write(format_str % (fraction, geneNum_median))
                # def format_float(x): return round(x / 100, 4)
                res_dict["fraction"].append(round(fraction, 1))
                res_dict["median_gene"].append(geneNum_median)

        return res_dict

    @utils.add_log
    def add_count_metrics(self, df):
        out_condi = " "
        out_condi += f"umi>={self.umi_cutoff} " if self.umi_cutoff > 0 else ""
        out_condi += f"read>={self.read_cutoff} " if self.read_cutoff > 0 else ""
        out_condi += f"gene>={self.gene_cutoff} " if self.gene_cutoff > 0 else ""
        self.add_help_content(
            name="Well:",
            content=f"The 1st column in the table is the well barcode sequence. Only output wells that meet the conditions({out_condi}).",
        )
        stats = df.describe()
        stats.columns = ["UMI", "Reads", "Genes"]
        for item in ["Reads", "UMI", "Genes"]:
            self.add_metric(
                name=f"Median {item} per Well",
                value=int(stats.loc["50%", item]),
                help_info="",
            )
        for item in ["Reads", "UMI", "Genes"]:
            self.add_metric(
                name=f"Mean {item} per Well",
                value=int(stats.loc["mean", item]),
                help_info="",
            )
        self.add_metric(
            name="Sequencing saturation",
            value=self.saturation,
            value_type="fraction",
            help_info="the fraction of read originating from an already-observed UMI.",
        )
        table_dict = self.get_table_dict(
            title="Detailed information per Well",
            table_id=self._table_id,
            df_table=df.reset_index(),
        )
        self.add_data(table_dict=table_dict)

    def add_plot_data(self):
        df_plot = pd.DataFrame(self.downsample_dict)
        df_plot = df_plot.rename(
            columns={"median_gene": "Median Genes", "fraction": "Reads Fraction"}
        )
        self.add_data(
            df_line=Line_plot(
                df_plot, x_title="Reads Fraction", y_title="Median Genes"
            ).get_plotly_div()
        )


def count(args):
    with Count(args, display_title="Wells") as runner:
        runner.run()


def get_opts_count(parser, sub_program):
    parser.add_argument(
        "--umi_cutoff",
        default=500,
        type=int,
        help="If the UMI number exceeds the threshold, it is considered a valid well and reported.",
    )
    parser.add_argument(
        "--gene_cutoff",
        default=0,
        type=int,
        help="If the gene number exceeds the threshold, it is considered a valid well and reported.",
    )
    parser.add_argument(
        "--read_cutoff",
        default=0,
        type=int,
        help="If the read number exceeds the threshold, it is considered a valid well and reported.",
    )
    super_count.get_opts_count(parser, sub_program)
