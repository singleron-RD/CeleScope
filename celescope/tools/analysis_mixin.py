import subprocess

import pandas as pd

import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH
from celescope.rna.mkref import parse_genomeDir_rna
from celescope.tools.step import Step


class AnalysisMixin(Step):
    """
    mixin class for analysis
    """

    def __init__(self, args, display_title=None):

        super().__init__(args, display_title=display_title)

        self.match_dir = None
        if hasattr(args, "match_dir") and args.match_dir:
            self.match_dir = args.match_dir
            self.read_match_dir()
        elif hasattr(args, "tsne_file") and args.tsne_file:
            tsne_df_file = args.tsne_file
            self.df_tsne = pd.read_csv(tsne_df_file, sep="\t")
            self.df_tsne.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            marker_df_file = args.marker_file
            self.read_format_df_marker(marker_df_file)
        else:
            self.match_dir = args.outdir + "/../"  # use self

    @utils.add_log
    def seurat(self, matrix_file, save_rds, genomeDir):
        app = ROOT_PATH + "/tools/run_analysis.R"
        genome = parse_genomeDir_rna(genomeDir)
        mt_gene_list = genome['mt_gene_list']
        cmd = (
            f'Rscript {app} '
            f'--sample {self.sample} '
            f'--outdir {self.outdir} '
            f'--matrix_file {matrix_file} '
            f'--mt_gene_list {mt_gene_list} '
            f'--save_rds {save_rds}'
        )
        self.seurat.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def auto_assign(self, type_marker_tsv):
        rds = f'{self.outdir}/{self.sample}.rds'
        app = ROOT_PATH + "/tools/auto_assign.R"
        cmd = (
            f'Rscript {app} '
            f'--rds {rds} '
            f'--type_marker_tsv {type_marker_tsv} '
            f'--outdir {self.outdir} '
            f'--sample {self.sample} '
        )
        self.auto_assign.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def read_format_df_marker(self):

        df_marker = pd.read_csv(self.df_marker_file, sep="\t")
        avg_logfc_col = "avg_log2FC"  # seurat 4
        if "avg_logFC" in df_marker.columns:  # seurat 2.3.4
            avg_logfc_col = "avg_logFC"
        df_marker = df_marker.loc[:,
                                  ["cluster", "gene", avg_logfc_col, "pct.1", "pct.2", "p_val_adj"]
                                  ]
        df_marker["cluster"] = df_marker["cluster"].apply(lambda x: f"cluster {x}")

        self.df_marker = df_marker

    def read_match_dir(self):
        """
        if match_dir is not self, should read match_dir at init
        if it is self, read at run_analysis - need to run seurat first
        """
        if self.match_dir:
            match_dict = utils.parse_match_dir(self.match_dir)
            tsne_df_file = match_dict['tsne_coord']
            self.df_tsne = pd.read_csv(tsne_df_file, sep="\t")
            self.df_tsne.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            self.df_marker_file = match_dict['markers']
            self.read_format_df_marker()
