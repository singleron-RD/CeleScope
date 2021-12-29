import subprocess

import pandas as pd

import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH
from celescope.rna.mkref import Mkref_rna
from celescope.tools.step import Step, s_common


def get_opts_analysis_mixin(parser, sub_program):
    
    parser.add_argument('--genomeDir', help='Required. Genome directory.', required=True)
    parser.add_argument('--save_rds', action='store_true', help='Write rds to disk.')
    parser.add_argument(
        '--type_marker_tsv',
        help="""A tsv file with header. If this parameter is provided, cell type will be annotated. Example:
```
cell_type	marker
Alveolar	"CLDN18,FOLR1,AQP4,PEBP4"
Endothelial	"CLDN5,FLT1,CDH5,RAMP2"
Epithelial	"CAPS,TMEM190,PIFO,SNTN"
Fibroblast	"COL1A1,DCN,COL1A2,C1R"
B_cell	"CD79A,IGKC,IGLC3,IGHG3"
Myeloid	"LYZ,MARCO,FCGR3A"
T_cell	"CD3D,TRBC1,TRBC2,TRAC"
LUAD	"NKX2-1,NAPSA,EPCAM"
LUSC	"TP63,KRT5,KRT6A,KRT6B,EPCAM"
```"""
    )
    if sub_program:
        parser.add_argument(
            '--matrix_file',
            help='Required. Matrix_10X directory from step count.',
            required=True,
        )
        # do not need all count diretory
        parser.add_argument("--tsne_file", help="match_dir t-SNE coord file.")
        parser.add_argument("--df_marker_file", help="match_dir df_marker_file.")
        parser = s_common(parser)

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
            self.df_tsne = pd.read_csv(tsne_df_file, sep="\t", index_col=0)
            self.df_tsne.index.rename("barcode", inplace=True)
            self.df_marker_file = args.df_marker_file
            self.read_format_df_marker()
        else:
            self.match_dir = args.outdir + "/../"  # use self

        # out
        self.rds = f'{self.out_prefix}.rds'
        self.mito_metric_file = f'{self.outdir}/mito.txt'

    @utils.add_log
    def seurat(self, matrix_file, save_rds, genomeDir):
        app = ROOT_PATH + "/tools/run_analysis.R"
        genome = Mkref_rna.parse_genomeDir(genomeDir)
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
    def add_mito_metric(self):
        with open(self.mito_metric_file, 'r') as f:
            for line in f:
                tmp = line.split(':')
                if tmp:
                    name, value = tmp[:2]
                    value = float(value)
                    self.add_metric(
                        name=name,
                        value=value,
                        display=f'{value}%'
                    )


    @utils.add_log
    def auto_assign(self, type_marker_tsv):
        app = ROOT_PATH + "/tools/auto_assign.R"
        cmd = (
            f'Rscript {app} '
            f'--rds {self.rds} '
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

    def run(self):
        pass
