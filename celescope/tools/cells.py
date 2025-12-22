import shutil
import os
import sys
import subprocess

import numpy as np
import pandas as pd

from celescope.tools.__init__ import (
    FILTERED_MATRIX_DIR_SUFFIX,
    RAW_MATRIX_DIR_SUFFIX,
    COUNTS_FILE_NAME,
)
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.tools.matrix import CountMatrix
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.rna.mkref import Mkref_rna


class Cells_metrics(Step):
    @utils.add_log
    def add_cells_metrics(
        self,
        n_cells,
        fraction_reads_in_cells,
        mean_used_reads_per_cell,
        median_umi_per_cell,
        total_genes,
        median_genes_per_cell,
        saturation,
        valid_reads,
    ):
        self.add_metric(
            name="Estimated Number of Cells",
            value=n_cells,
            help_info="The number of barcodes considered as cell-associated.",
        )

        self.add_metric(
            name="Fraction Reads in Cells",
            value=fraction_reads_in_cells,
            value_type="fraction",
            help_info="the fraction of uniquely-mapped-to-transcriptome reads with cell-associated barcodes",
        )

        mean_reads_per_cell = valid_reads // n_cells
        self.add_metric(
            name="Mean Reads per Cell",
            value=mean_reads_per_cell,
            help_info="the number of Valid Reads divided by Estimated Number of Cells",
        )

        self.add_metric(
            name="Mean Used Reads per Cell",
            value=mean_used_reads_per_cell,
            help_info="the number of uniquely-mapped-to-transcriptome reads per cell-associated barcode",
        )

        self.add_metric(
            name="Median UMI per Cell",
            value=median_umi_per_cell,
            help_info="the median number of UMI counts per cell-associated barcode",
        )

        self.add_metric(
            name="Total Genes",
            value=total_genes,
            help_info="the number of genes with at least one UMI count in any cell",
        )

        self.add_metric(
            name="Median Genes per Cell",
            value=median_genes_per_cell,
            help_info="the median number of genes detected per cell-associated barcode",
        )

        self.add_metric(
            name="Saturation",
            value=saturation,
            value_type="fraction",
            help_info="the fraction of read originating from an already-observed UMI. ",
        )

    def run(self):
        pass


class Cells(Cells_metrics):
    def __init__(self, args, display_title=None):
        os.chdir(args.root_dir)
        args.outdir = f"{args.sample}/cells"
        super().__init__(args, display_title=display_title)
        self.raw_matrix = f"{self.outs_dir}/{RAW_MATRIX_DIR_SUFFIX}"
        self.old_filtered_matrix = f"{self.outs_dir}/{FILTERED_MATRIX_DIR_SUFFIX}"
        self.counts_file = f"{self.outs_dir}/{COUNTS_FILE_NAME}"

        # out
        self.filter_matrix = f"{self.outdir}/{FILTERED_MATRIX_DIR_SUFFIX}"
        self.default_filter_matrix = (
            f"{self.outs_dir}/default_{FILTERED_MATRIX_DIR_SUFFIX}"
        )
        self.default_report_html = f"{self.outdir}/../default_{self.sample}_report.html"

        if os.path.exists(self.default_filter_matrix):
            self.old_filtered_matrix = self.default_filter_matrix

        self.outs = [self.filter_matrix]

    @utils.add_log
    def force_cells(self):
        raw = CountMatrix.from_matrix_dir(self.raw_matrix)
        df_counts = pd.read_csv(self.counts_file, index_col=0, header=0, sep="\t")
        bcs = list(df_counts.head(self.args.force_cells).index)
        filtered = raw.slice_matrix_bc(bcs)
        self.add_metric(
            name="Force cells",
            value=self.args.force_cells,
            show=False,
        )
        return filtered

    @utils.add_log
    def filter_min_gene(self, filtered: CountMatrix) -> CountMatrix:
        bc_geneNum, _ = filtered.get_bc_geneNum()
        bc_indices = [bc for bc in bc_geneNum if bc_geneNum[bc] >= self.args.min_gene]
        n_filtered_cells = len(bc_geneNum) - len(bc_indices)
        self.add_metric(
            name="Minimum gene number in cells",
            value=self.args.min_gene,
            show=False,
        )
        self.add_metric(
            name="Number of cells filtered by min_gene",
            value=n_filtered_cells,
            show=False,
        )
        return filtered.slice_matrix_bc(bc_indices)

    @utils.add_log
    def metrics_report(self, filtered: CountMatrix):
        p_default_dict = {
            self.old_filtered_matrix: self.default_filter_matrix,
            self.report_html: self.default_report_html,
        }
        for p, default_p in p_default_dict.items():
            if not os.path.exists(default_p):
                shutil.move(p, default_p)

        self.add_slot_step(
            slot="metrics",
            step_name="default_cells",
            val=self.old_step_dict["metrics"],
        )

        self.add_slot_step(
            slot="data",
            step_name="default_cells",
            val=self.old_step_dict["data"],
        )
        bcs = filtered.get_barcodes()
        n_cells = len(bcs)

        df_counts = pd.read_csv(self.counts_file, index_col=0, header=0, sep="\t")
        reads_total = df_counts["countedU"].sum()
        reads_cell = df_counts.loc[bcs, "countedU"].sum()
        fraction_reads_in_cells = float(reads_cell / reads_total)
        mean_used_reads_per_cell = int(reads_cell // len(bcs))
        median_umi_per_cell = int(df_counts.loc[bcs, "UMI"].median())

        bc_geneNum, total_genes = filtered.get_bc_geneNum()
        median_genes_per_cell = int(np.median(list(bc_geneNum.values())))

        saturation = self.old_step_dict["metrics"]["Saturation"] / 100
        df_counts.loc[:, "mark"] = "UB"
        df_counts.loc[bcs, "mark"] = "CB"
        df_counts.to_csv(self.counts_file, sep="\t", index=True)
        valid_reads = (
            self.old_step_dict["metrics"]["Mean Reads per Cell"]
            * self.old_step_dict["metrics"]["Estimated Number of Cells"]
        )
        self.add_cells_metrics(
            n_cells,
            fraction_reads_in_cells,
            mean_used_reads_per_cell,
            median_umi_per_cell,
            total_genes,
            median_genes_per_cell,
            saturation,
            valid_reads,
        )
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.counts_file))

    @utils.add_log
    def filter_mito(self, filtered: CountMatrix) -> CountMatrix:
        if not self.args.genomeDir:
            raise ValueError("--genomeDir must be provided when --max_mito is set.")
        mt_gene_list = Mkref_rna.get_config(self.args.genomeDir)["files"][
            "mt_gene_list"
        ]
        if not mt_gene_list:
            raise FileNotFoundError(
                f"mt_gene_list not exist in genomeDir: {self.args.genomeDir}"
            )
        mito_genes, _ = utils.read_one_col(mt_gene_list)
        bc_fraction = filtered.get_genes_fraction(mito_genes)
        bc_indices = (bc_fraction <= self.args.max_mito).nonzero()[1]
        n_filtered_cells = len(filtered.get_barcodes()) - len(bc_indices)
        filtered = filtered.slice_matrix(bc_indices)
        self.add_metric(
            name="Maximum fraction of mitocondrial UMI in cells",
            value=self.args.max_mito,
            show=False,
        )
        self.add_metric(
            name="Number of cells filtered by max_mito",
            value=n_filtered_cells,
            show=False,
        )
        return filtered

    @utils.add_log
    def soloCellFilter(self, filtered: CountMatrix) -> CountMatrix:
        temp_raw = f"{self.outdir}/raw"
        cmd = f"cp -r {self.raw_matrix} {self.outdir}; gunzip {temp_raw}/*.gz"  # STAR do not recognize gzip matrix
        if not os.path.exists(temp_raw):
            subprocess.check_call(cmd, shell=True)
        temp_filtered = f"./{self.outdir}/filtered/"
        cmd = f"STAR --runMode soloCellFiltering {temp_raw} {temp_filtered} --soloCellFilter {self.args.soloCellFilter}"
        subprocess.check_call(cmd, shell=True)
        cmd = f"gzip {temp_filtered}/*"
        subprocess.check_call(cmd, shell=True)
        filtered = CountMatrix.from_matrix_dir(temp_filtered)
        self.add_metric(
            name="soloCellFilter",
            value=self.args.soloCellFilter,
            show=False,
        )
        return filtered

    @utils.add_log
    def force_barcode(self):
        raw = CountMatrix.from_matrix_dir(self.raw_matrix)
        bcs = utils.one_col_to_list(self.args.barcode)
        filtered = raw.slice_matrix_bc(bcs)
        return filtered

    @utils.add_log
    def run(self):
        if self.args.max_mito > 1.0:
            sys.exit("max_mito should be less than 1.0")
        filtered = CountMatrix.from_matrix_dir(self.old_filtered_matrix)
        if self.args.barcode:
            filtered = self.force_barcode()
        elif self.args.force_cells > 0:
            filtered = self.force_cells()
        elif self.args.soloCellFilter:
            filtered = self.soloCellFilter(filtered)
        if self.args.max_mito < 1.0:
            filtered = self.filter_mito(filtered)
        if self.args.min_gene > 0:
            filtered = self.filter_min_gene(filtered)

        filtered.to_matrix_dir(self.filter_matrix)
        self.metrics_report(filtered)


def cells(args):
    with Cells(args) as runner:
        runner.run()


def get_opts_cells(parser, sub_program=True):
    if sub_program:
        parser.add_argument(
            "--force_cells",
            help="Force to use this number of cells.",
            default=0,
            type=int,
        )
        parser.add_argument(
            "--soloCellFilter",
            help="The same as the argument in STARsolo. Ignored when --force_cells is set.",
            default="",
        )
        parser.add_argument(
            "--root_dir",
            help="The root directory of CeleScope runs.",
            default="./",
        )
        parser.add_argument(
            "--max_mito",
            help="Maximum mitocondrial fraction in a cell.",
            default="1.0",
            type=float,
        )
        parser.add_argument(
            "--min_gene",
            help="Minimum gene number in a cell.",
            default=0,
            type=int,
        )
        parser.add_argument(
            "--genomeDir",
            help=HELP_DICT["genomeDir"],
        )
        parser.add_argument(
            "--barcode",
            help="User defined barcode file. One barcode per line.",
        )
        parser = s_common(parser)
