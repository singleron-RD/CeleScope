from celescope.tools.matrix import CountMatrix
from celescope.tools.starsolo import (
    Starsolo as tools_Starsolo,
    get_opts_starsolo as tools_opts,
    Mapping,
    Demultiplexing,
)
from celescope.space.utils import Spatial
from celescope.tools.__init__ import COUNTS_FILE_NAME
import pandas as pd
from celescope.tools.emptydrop_cr import get_plot_elements
import numpy as np
from celescope.tools.utils import add_log
from celescope.tools.step import Step


class Starsolo(tools_Starsolo):
    def __init__(self, args):
        super().__init__(args)
        if self.chemistry == "space-ffpe":
            self.extra_starsolo_args += " --soloStrand Reverse --clip5pNbases 44 "
            args.fq2 = args.fq1  # single end

    @add_log
    def keep_barcodes(self):
        matrix = CountMatrix.from_matrix_dir(self.raw_matrix)
        in_tissue_barcodes = Spatial(self.args.spatial).get_in_tissue_barcodes()
        filtered = matrix.slice_matrix_bc(in_tissue_barcodes)
        filtered.to_matrix_dir(self.filtered_matrix)
        return filtered

    def run(self):
        self.run_starsolo()
        filtered = self.keep_barcodes()
        self.gzip_matrix()
        q30_cb, q30_umi = self.get_Q30_cb_UMI()
        return q30_cb, q30_umi, filtered


class Spots(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        solo_dir = (
            f"{self.outdir}/{self.sample}_Solo.out/{self.args.report_soloFeature}"
        )
        self.summary_file = f"{solo_dir}/Summary.csv"
        self.counts_file = f"{self.outs_dir}/{COUNTS_FILE_NAME}"

    @add_log
    def parse_summary(self):
        df = pd.read_csv(self.summary_file, index_col=0, header=None)
        s = df.iloc[:, 0]
        saturation = float(s["Sequencing Saturation"])
        n_reads = int(s["Number of Reads"])
        q30_RNA = float(s["Q30 Bases in RNA read"])

        return n_reads, q30_RNA, saturation

    def run(
        self,
        filtered: CountMatrix,
    ):
        df_counts = pd.read_csv(self.counts_file, index_col=0, header=0, sep="\t")
        df_counts = df_counts.loc[df_counts["countedU"] > 0]
        reads_total = df_counts["countedU"].sum()
        bcs = set(df_counts.index).intersection(filtered.get_barcodes())
        n_spots = len(bcs)
        reads_spot = df_counts.loc[bcs, "countedU"].sum()
        fraction_reads_in_spots = float(reads_spot / reads_total)
        mean_used_reads_per_spot = int(reads_spot // len(bcs))
        median_umi_per_spot = int(df_counts.loc[bcs, "UMI"].median())

        bc_geneNum, total_genes = filtered.get_bc_geneNum()
        median_genes_per_spot = int(np.median(list(bc_geneNum.values())))

        df_counts.loc[:, "mark"] = "UB"
        df_counts.loc[bcs, "mark"] = "CB"
        df_counts.fillna(0, inplace=True)
        df_counts = df_counts.astype({"UMI": int, "countedU": int})
        df_counts.to_csv(self.counts_file, sep="\t", index=True)

        self.add_metric(
            "In Tissue Spots",
            n_spots,
            help_info="Number of spots in tissue determined based on the image",
        )
        self.add_metric(
            "Fraction of Reads in Spots",
            fraction_reads_in_spots,
            value_type="fraction",
            help_info="Fraction of reads which were mapped to a in tissue spot",
        )
        self.add_metric(
            "Mean Used Reads per Spot",
            mean_used_reads_per_spot,
            help_info="The number of uniquely-mapped-to-transcriptome reads per in tissue spot",
        )
        self.add_metric(
            "Median UMI per Spot",
            median_umi_per_spot,
            help_info="Median UMI count per spot",
        )
        self.add_metric(
            "Median Genes per Spot",
            median_genes_per_spot,
            help_info="Median number of genes per spot",
        )
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.counts_file))

        n_reads, q30_RNA, saturation = self.parse_summary()
        self.add_metric(
            "Saturation",
            saturation,
            value_type="fraction",
            help_info="the fraction of read originating from an already-observed UMI.",
        )

        return n_reads, q30_RNA


def starsolo(args):
    with Starsolo(args) as runner:
        q30_cb, q30_umi, filtered = runner.run()

    with Mapping(args) as runner:
        valid_reads, corrected = runner.run()

    with Spots(args) as runner:
        n_reads, q30_RNA = runner.run(filtered)

    with Demultiplexing(args) as runner:
        runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)


def get_opts_starsolo(parser, sub_program=True):
    tools_opts(parser, sub_program)
    parser.set_defaults(
        outFilterMatchNmin=30,
        soloCellFilter="None",
    )
    if sub_program:
        parser.add_argument(
            "--spatial",
            help="spatial directory.",
            required=True,
        )
    return parser
