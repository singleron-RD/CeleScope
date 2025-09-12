import subprocess
import sys

import numpy as np
import pandas as pd

from celescope.tools.__init__ import COUNTS_FILE_NAME
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools.starsolo import (
    Demultiplexing,
    Mapping as ToolsMapping,
    create_solo_str,
    create_soloFeatures,
    get_opts_starsolo as tools_opts,
)
from celescope.tools.starsolo import (
    Starsolo as tools_Starsolo,
)
from celescope.tools.matrix import CountMatrix
from celescope.tools import utils
from celescope.tools.step import Step

SAM_attributes = "NH HI nM AS CR UR CB UB GX GN "


def get_well_barcode(
    bc_file: str,
) -> dict[int, str]:
    barcodes = utils.one_col_to_list(bc_file)
    return {i: x for i, x in enumerate(barcodes, start=1)}


def get_barcode_sample(
    bc_file: str,
    well_sample_file: str,
) -> dict[str, str]:
    well_barcode = get_well_barcode(bc_file)
    well_sample = utils.two_col_to_dict(well_sample_file)
    barcode_sample = {
        well_barcode[well]: sample for well, sample in well_sample.items()
    }
    return barcode_sample


class Starsolo(tools_Starsolo):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        self.well_barcode = get_well_barcode(self.bc[0])
        self.barcode_sample = get_barcode_sample(self.bc[0], args.well_sample)

        self.tsv_matrix_file = f"{self.out_prefix}_matrix.tsv.gz"
        self.outs.append(self.tsv_matrix_file)

    def run_starsolo(self):
        if "--outSAMunmapped Within" not in self.extra_starsolo_args:
            self.extra_starsolo_args += " --outSAMunmapped Within "
        soloFeatures = create_soloFeatures(
            self.args.soloFeatures, self.args.report_soloFeature
        )
        cmd = create_solo_str(
            pattern_args=self.pattern_args,
            whitelist_args=self.whitelist_args,
            outFileNamePrefix=self.out_prefix + "_",
            fq1=self.args.fq1,
            fq2=self.args.fq2,
            genomeDir=self.args.genomeDir,
            soloCellFilter="None",  # only use raw matrix
            runThreadN=self.args.thread,
            clip3pAdapterSeq=self.args.adapter_3p,
            outFilterMatchNmin=self.args.outFilterMatchNmin,
            soloFeatures=soloFeatures,
            outSAMtype=self.args.outSAMtype,
            outSAMattributes=self.outSAMattributes,
            soloCBmatchWLtype=self.args.soloCBmatchWLtype,
            limitBAMsortRAM=self.args.limitBAMsortRAM,
            extra_starsolo_args=self.extra_starsolo_args,
        )
        if self.chemistry == "bulk_rna-bulk_vdj_match":
            cmd += "--soloStrand Reverse \\\n"
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)
        cmd = f"chmod -R 755 {self.solo_out_dir}"
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def keep_barcodes(self):
        """take in raw matrix, only keep barcodes in the input file, and convert barcodes to sample names"""
        matrix = CountMatrix.from_matrix_dir(self.raw_matrix)

        filtered = matrix.slice_matrix_bc(self.barcode_sample.keys())
        filtered.to_matrix_dir(self.filtered_matrix)
        samples = [self.barcode_sample[bc] for bc in filtered.get_barcodes()]
        converted = CountMatrix(filtered.get_features(), samples, filtered.get_matrix())
        df = converted.to_df()
        df.to_csv(self.tsv_matrix_file, sep="\t")
        return filtered

    def run(self):
        self.run_starsolo()
        filtered = self.keep_barcodes()
        self.gzip_matrix()
        q30_cb, q30_umi = self.get_Q30_cb_UMI()
        return q30_cb, q30_umi, filtered, self.barcode_sample, self.well_barcode


class Mapping(ToolsMapping):
    # only add metrics
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    @utils.add_log
    def parse_cellReadsStats(self) -> tuple[dict[str, int], pd.DataFrame]:
        df = pd.read_csv(self.cellReadsStats, sep="\t", header=0, index_col=0)
        df = df.iloc[1:,]  # skip first line cb not pass whitelist
        df_count = df.loc[:, ["nUMIunique", "countedU"]]  # keep dataframe format
        df_count.rename(columns={"nUMIunique": "UMI"}, inplace=True)
        df_count.sort_values(by="UMI", ascending=False, inplace=True)
        cbs = CountMatrix.read_barcodes(self.filtered_matrix)
        df_count["mark"] = "UB"
        for cb in cbs:
            df_count.loc[cb, "mark"] = "CB"
        df_count.to_csv(self.counts_file, sep="\t", index=True)

        df = df.loc[
            :,
            [
                "cbMatch",
                "cbPerfect",
                "genomeU",
                "genomeM",
                "exonic",
                "intronic",
                "exonicAS",
                "intronicAS",
                "countedU",
            ],
        ]
        metrics = df.sum().to_dict()
        df = df.loc[
            :,
            [
                "cbMatch",
                "genomeU",
                "genomeM",
            ],
        ]
        df.rename(
            columns={
                "cbMatch": "Reads",
                "genomeU": "Unique-mapping Reads",
                "genomeM": "Multiple-mapping Reads",
            },
            inplace=True,
        )
        df["Unique-mapping Reads"] = np.where(
            df["Reads"] != 0,
            (df["Unique-mapping Reads"] / df["Reads"] * 100).round(2).astype(str) + "%",
            np.nan,
        )
        df["Multiple-mapping Reads"] = np.where(
            df["Reads"] != 0,
            (df["Multiple-mapping Reads"] / df["Reads"] * 100).round(2).astype(str)
            + "%",
            np.nan,
        )
        return metrics, df

    @utils.add_log
    def run(self):
        metrics, df = self.parse_cellReadsStats()
        return self.add_metrics_to_report(metrics), df


class Cells(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        solo_dir = (
            f"{self.outdir}/{self.sample}_Solo.out/{self.args.report_soloFeature}"
        )
        self.summary_file = f"{solo_dir}/Summary.csv"
        self.counts_file = f"{self.outs_dir}/{COUNTS_FILE_NAME}"

    @utils.add_log
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
        barcode_sample: dict,
        well_barcode: dict,
        df_metrics: pd.DataFrame,
    ):
        df_counts = pd.read_csv(self.counts_file, index_col=0, header=0, sep="\t")
        reads_total = df_counts["countedU"].sum()
        bcs = filtered.get_barcodes()
        n_cells = len(bcs)
        reads_cell = df_counts.loc[bcs, "countedU"].sum()
        fraction_reads_in_cells = float(reads_cell / reads_total)
        mean_used_reads_per_cell = int(reads_cell // len(bcs))
        median_umi_per_cell = int(df_counts.loc[bcs, "UMI"].median())

        bc_geneNum, total_genes = filtered.get_bc_geneNum()
        median_genes_per_cell = int(np.median(list(bc_geneNum.values())))

        df_counts.loc[:, "mark"] = "UB"
        df_counts.loc[bcs, "mark"] = "CB"
        df_counts.fillna(0, inplace=True)
        df_counts = df_counts.astype({"UMI": int, "countedU": int})
        df_counts.to_csv(self.counts_file, sep="\t", index=True)

        self.add_metric(
            "Number of Wells",
            n_cells,
            help_info="number of barcodes with at least one UMI",
        )
        self.add_metric(
            "Fraction of Reads in Wells",
            fraction_reads_in_cells,
            value_type="fraction",
            help_info="Fraction of reads which were mapped to a barcode",
        )
        self.add_metric(
            "Mean Used Reads per Well",
            mean_used_reads_per_cell,
            help_info="The number of uniquely-mapped-to-transcriptome reads per well",
        )
        self.add_metric(
            "Median UMI per Well",
            median_umi_per_cell,
            help_info="Median UMI count per barcode",
        )
        self.add_metric(
            "Median Genes per Well",
            median_genes_per_cell,
            help_info="Median number of genes per barcode",
        )
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.counts_file))

        n_reads, q30_RNA, saturation = self.parse_summary()
        self.add_metric(
            "Saturation",
            saturation,
            value_type="fraction",
            help_info="the fraction of read originating from an already-observed UMI.",
        )

        # table
        df_cells = df_counts[df_counts["mark"] == "CB"]
        df_cells["Sample"] = df_cells.index.map(lambda x: barcode_sample[x])
        df_cells["Genes"] = df_cells.index.map(lambda x: bc_geneNum.get(x, 0))
        df_cells["Barcode"] = df_cells.index
        barcode_well = {v: k for k, v in well_barcode.items()}
        df_cells["Well"] = df_cells["Barcode"].map(lambda x: barcode_well[x])
        df_cells = df_cells.merge(
            df_metrics, left_on="Barcode", right_index=True, how="left"
        )
        df_cells = df_cells.loc[
            :,
            [
                "Well",
                "Sample",
                "Barcode",
                "Reads",
                "Unique-mapping Reads",
                "Multiple-mapping Reads",
                "UMI",
                "Genes",
            ],
        ]

        table_dict = self.get_table_dict(
            title="Metrics for a Single Well",
            table_id=self.assay,
            df_table=df_cells,
        )
        self.add_data(table_dict=table_dict)

        return n_reads, q30_RNA


def starsolo(args):
    with Starsolo(args) as runner:
        q30_cb, q30_umi, filtered, barcode_sample, well_barcode = runner.run()

    with Mapping(args) as runner:
        (valid_reads, corrected), df_metrics = runner.run()

    with Cells(args, display_title="Wells") as runner:
        n_reads, q30_RNA = runner.run(
            filtered, barcode_sample, well_barcode, df_metrics
        )

    with Demultiplexing(args) as runner:
        runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)


def get_opts_starsolo(parser, sub_program=True):
    tools_opts(parser, sub_program)
    parser.add_argument(
        "--well_sample",
        help="tsv file of well numbers and sample names. The first column is well numbers and the second column is sample names.",
        required=True,
    )
    return parser
