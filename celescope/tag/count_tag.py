"""
assign cell identity based on SNR and UMI_min
"""

from celescope.__init__ import ROOT_PATH
from celescope.tools.step import Step, s_common
import celescope.tools.utils as utils
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import matplotlib

matplotlib.use('Agg')


def get_opts_count_tag(parser, sub_program):
    parser.add_argument(
        "--UMI_min",
        help="Default='auto'. Minimum UMI threshold. Cell barcodes with valid UMI < UMI_min are classified as *undeterminded*.",
        default="auto"
    )
    parser.add_argument(
        "--dim",
        help="Default=1. Tag dimentions. Usually we use 1-dimentional tag.",
        default=1
    )
    parser.add_argument(
        "--SNR_min",
        help="""Default='auto'. Minimum signal-to-noise ratio. 
Cell barcodes with UMI >=UMI_min and SNR < SNR_min are classified as *multiplet*. """,
        default="auto"
    )
    parser.add_argument("--combine_cluster",
                        help="Conbine cluster tsv file.", default=None)
    parser.add_argument(
        "--coefficient",
        help="""Default=0.1. If `SNR_min` is 'auto', minimum signal-to-noise ratio is calulated as 
`SNR_min = max(median(SNRs) * coefficient, 2)`. 
Smaller `coefficient` will cause less *multiplet* in the tag assignment.""",
        default=0.1
    )
    if sub_program:
        parser.add_argument("--read_count_file", help="Tag read count file.", required=True)
        parser.add_argument("--match_dir", help="Match celescope scRNA-Seq directory.")
        parser.add_argument("--matrix_dir", help="Match celescope scRNA-Seq matrix directory.")
        parser.add_argument("--tsne_file", help="t-SNE coord file.")

        s_common(parser)


def count_tag(args):

    with Count_tag(args, display_title="Cells") as runner:
        runner.run()


class Count_tag(Step):
    """
    Features
    - Assign tag to each cell barcode and summarize.

    Output

    - `{sample}_umi_tag.tsv` 

        `first column` cell barcode  
        `last column`  assigned tag  
        `columns between first and last` UMI count for each tag 

    - `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation

    - `{sample}_cluster_count.tsv` cell barcode number assigned to *undeterminded*, *multiplet* and *each tag*

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.read_count_file = args.read_count_file
        self.UMI_min = args.UMI_min
        self.SNR_min = args.SNR_min
        self.combine_cluster = args.combine_cluster
        self.dim = int(args.dim)
        self.coefficient = float(args.coefficient)

        # read
        self.df_read_count = pd.read_csv(self.read_count_file, sep="\t", index_col=0)

        if args.match_dir:
            match_dict = utils.parse_match_dir(args.match_dir)
            self.match_barcode = match_dict['match_barcode']
            self.n_match_barcode = match_dict['n_match_barcode']
            self.tsne_file = match_dict['tsne_coord']
            self.matrix_dir = match_dict['matrix_dir']
        elif args.matrix_dir:
            df_barcode = pd.read_csv(f'{args.matrix_dir}/barcodes.tsv', header=None)
            self.match_barcode = df_barcode[0].tolist()
            self.n_match_barcode = len(self.match_barcode)
            self.tsne_file = args.tsne_file
            self.matrix_dir = args.matrix_dir
        else:
            raise ValueError("--match_dir or --matrix_dir is required.")

        # init
        self.no_noise = False

        # out files
        self.UMI_tag_file = f'{self.outdir}/{self.sample}_umi_tag.tsv'
        self.tsne_tag_file = f'{self.outdir}/{self.sample}_tsne_tag.tsv'
        self.cluster_count_file = f'{self.outdir}/{self.sample}_cluster_count.tsv'
        self.cluster_plot = f'{self.outdir}/{self.sample}_cluster_plot.pdf'
        if self.combine_cluster:
            self.combine_cluster_count_file = f'{self.outdir}/{self.sample}_combine_cluster_count.tsv'
            self.combine_cluster_plot = f'{self.outdir}/{self.sample}_combine_cluster_plot.pdf'

    @staticmethod
    def get_UMI(row):
        return row.sum()

    @staticmethod
    def get_UMI_min(df_cell_UMI, UMI_min):
        if UMI_min == "auto":
            UMI_min1 = np.percentile(df_cell_UMI.sum(axis=1), 5)
            UMI_min2 = np.median(df_cell_UMI.sum(axis=1)) / 10
            UMI_min = int(min(UMI_min1, UMI_min2))
            UMI_min = max(UMI_min, 1)
            return UMI_min
        else:
            return int(UMI_min)

    @staticmethod
    def get_SNR(row, dim):
        row_sorted = sorted(row, reverse=True)
        noise = row_sorted[dim]
        signal = row_sorted[dim - 1]
        if signal == 0:
            return 0
        if noise == 0:
            return np.inf
        return float(signal) / noise

    @utils.add_log
    def get_SNR_min(self, df_cell_UMI, SNR_min, UMI_min):
        UMIs = df_cell_UMI.apply(Count_tag.get_UMI, axis=1)
        df_valid_cell_UMI = df_cell_UMI[UMIs >= UMI_min]
        if SNR_min == "auto":
            # no noise
            if df_valid_cell_UMI.shape[1] <= self.dim:
                Count_tag.get_SNR_min.logger.warning('*** No NOISE FOUND! ***')
                self.no_noise = True
                return 0
            SNRs = df_valid_cell_UMI.apply(Count_tag.get_SNR, dim=self.dim, axis=1)
            if np.median(SNRs) == np.inf:
                return 10
            return max(np.median(SNRs) * self.coefficient, 2)
        else:
            return float(SNR_min)

    @staticmethod
    def tag_type(row, UMI_min, SNR_min, dim, no_noise=False):
        if no_noise:
            SNR = 1
        else:
            SNR = Count_tag.get_SNR(row, dim)
        UMI = Count_tag.get_UMI(row)
        if UMI < UMI_min:
            return "Undetermined"
        if SNR < SNR_min:
            return "Multiplet"
        # get tag
        signal_tags = sorted(row.sort_values(ascending=False).index[0:dim])
        signal_tags_str = "_".join(signal_tags)
        return signal_tags_str

    @utils.add_log
    def write_and_plot(self, df, column_name, count_file, plot_file):
        df_count = df.groupby(["tag", column_name]).size().unstack()
        df_count.fillna(0, inplace=True)
        df_count.to_csv(count_file, sep="\t")
        df_percent = df_count / df_count.sum()
        df_plot = df_percent.stack().reset_index()
        df_plot.rename({0: "percent"}, axis=1, inplace=True)

        # plot
        colors = list(matplotlib.colors.cnames.keys())
        fig, ax = plt.subplots(figsize=(20, 10))
        types = df_plot["tag"].drop_duplicates()
        margin_bottom = np.zeros(len(df_plot[column_name].drop_duplicates()))

        for num, tag_type in enumerate(types):
            try:
                color = colors[num * 3 + 1]
            except IndexError:
                print("not enough colors")
                return
            values = list(df_plot.loc[df_plot["tag"] == tag_type, "percent"])
            df_plot[df_plot['tag'] == tag_type].plot.bar(
                x=column_name, y='percent', ax=ax, stacked=True,
                bottom=margin_bottom, label=tag_type, color=color)

            margin_bottom += values
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title("tag fraction")
        fig.savefig(plot_file)

    @utils.add_log
    def run(self):

        mapped_read = int(self.df_read_count['read_count'].sum())

        # in cell
        df_read_count_in_cell = self.df_read_count[self.df_read_count.index.isin(self.match_barcode)]
        mapped_read_in_cell = int(df_read_count_in_cell['read_count'].sum())
        self.add_metric(
            name='Mapped Reads in Cells',
            value=mapped_read_in_cell,
            total=mapped_read,
            help_info="Mapped reads with scRNA-Seq cell barcode"
        )

        # UMI
        tag_name = df_read_count_in_cell.columns[0]
        df_UMI_in_cell = df_read_count_in_cell.reset_index().groupby([
            'barcode', tag_name]).agg({'UMI': 'count'})
        df_UMI_in_cell = df_UMI_in_cell.reset_index()
        df_UMI_in_cell = df_UMI_in_cell.pivot(
            index='barcode', columns=tag_name, values='UMI')
        df_cell = pd.DataFrame(index=self.match_barcode)
        df_UMI_cell = pd.merge(
            df_cell,
            df_UMI_in_cell,
            how="left",
            left_index=True,
            right_index=True
        )

        # fillna
        df_UMI_cell.fillna(0, inplace=True)
        df_UMI_cell = df_UMI_cell.astype(int)

        # UMI
        UMIs = df_UMI_cell.apply(sum, axis=1)
        umi_median = round(np.median(UMIs), 2)
        umi_mean = round(np.mean(UMIs), 2)
        self.add_metric(
            name='Median UMI per Cell',
            value=umi_median,
            help_info="Median UMI per scRNA-Seq cell barcode"
        )

        self.add_metric(
            name='Mean UMI per Cell',
            value=umi_mean,
            help_info="Mean UMI per scRNA-Seq cell barcode"
        )

        UMI_min = Count_tag.get_UMI_min(df_UMI_cell, self.UMI_min)
        Count_tag.run.logger.info(f'UMI_min: {UMI_min}')
        SNR_min = self.get_SNR_min(df_UMI_cell, self.SNR_min, UMI_min)
        Count_tag.run.logger.info(f'SNR_min: {SNR_min}')
        df_UMI_cell["tag"] = df_UMI_cell.apply(
            Count_tag.tag_type, UMI_min=UMI_min, SNR_min=SNR_min, dim=self.dim, no_noise=self.no_noise, axis=1)
        df_UMI_cell.to_csv(self.UMI_tag_file, sep="\t")

        df_tsne = pd.read_csv(self.tsne_file, sep="\t", index_col=0)
        df_tsne_tag = pd.merge(
            df_tsne,
            df_UMI_cell,
            how="left",
            left_index=True,
            right_index=True)

        if self.combine_cluster:
            df_combine_cluster = pd.read_csv(
                self.combine_cluster, sep="\t", header=None)
            df_combine_cluster.columns = ["cluster", "combine_cluster"]
            df_tsne_combine_cluster_tag = pd.merge(
                df_tsne_tag, df_combine_cluster,
                on=["cluster"], how="left", left_index=True).set_index(df_tsne_tag.index)
            df_tsne_combine_cluster_tag.to_csv(self.tsne_tag_file, sep="\t")
        else:
            df_tsne_tag.to_csv(self.tsne_tag_file, sep="\t")

        self.write_and_plot(
            df=df_tsne_tag,
            column_name="cluster",
            count_file=self.cluster_count_file,
            plot_file=self.cluster_plot,
        )

        if self.combine_cluster:
            self.write_and_plot(
                df=df_tsne_combine_cluster_tag,
                column_name="combine_cluster",
                count_file=self.combine_cluster_count_file,
                plot_file=self.combine_cluster_plot
            )

        sr_tag_count = df_UMI_cell["tag"].value_counts()  # series(index:tag name, value:tag count)
        for tag_name in ("Undetermined", "Multiplet"):
            if tag_name in sr_tag_count:
                self.add_metric(
                    name=tag_name + ' Cells',
                    value=int(sr_tag_count[tag_name]),
                    total=self.n_match_barcode,
                )
                sr_tag_count.drop(tag_name, inplace=True)
        for tag_name in sorted(sr_tag_count.index):
            self.add_metric(
                name=tag_name + ' Cells',
                value=int(sr_tag_count[tag_name]),
                total=self.n_match_barcode,
            )

        # seurat hashtag
        if self.debug:
            self.seurat_hashtag()

    @utils.add_log
    def seurat_hashtag(self):
        app = f'{ROOT_PATH}/tag/seurat_hashtag.R'
        cmd = (
            f'Rscript {app} '
            f'--outdir {self.outdir} '
            f'--sample {self.sample} '
            f'--umi_tag {self.UMI_tag_file} '
            f'--matrix_10X {self.matrix_dir} '
        )
        self.debug_subprocess_call(cmd)
