import pandas as pd
import glob
import subprocess
import os
import celescope
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import warnings
from celescope.tools import utils
from celescope.tools.step import s_common, Step
matplotlib.use('Agg')
warnings.filterwarnings("ignore")

TOOLS_DIR = os.path.dirname(celescope.tools.__file__)

CELL_TYPE_DICT = {
    'TCR':['T_cells', 'NKT_cells', 'T cells', 'NK T cells', 'Tcells'],
    'BCR':['Plasma_cells', 'B_cells', 'Mature_B_cell', 'Plasma cells', 'B cells', 'Bcells'],
    }


class Mapping(Step):
    """
    ## Features

    - Output TSNE-plot of Assembled T/B Cells.

    ## Output
    - `06.mapping/{sample}_mapping.pdf` TSNE-plot of Assembled Cells.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        
        self.seqtype = args.seqtype
        self.match_dir = args.match_dir
        if self.match_dir != 'None':
            self.contig_file = glob.glob(f'{args.match_out}/matched_contig_annotations.csv')[0]
        self.celltype_set = CELL_TYPE_DICT[self.seqtype]

        self.rds = glob.glob(f'{self.match_dir}/06.analysis/*.rds')
        self.assign_file = glob.glob(f'{self.match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')
        self.h5ad = glob.glob(f'{self.match_dir}/06.analysis/*.h5ad')


    @staticmethod
    @utils.add_log
    def run_mapping(rds, contig, sample, outdir, assign):
        cmd = (
            f'Rscript {TOOLS_DIR}/VDJmapping.R '
            f'--rds {rds} '
            f'--VDJ {contig} '
            f'--sample {sample} '
            f'--outdir {outdir} '
            f'--assign_file {assign} '
            '2>&1 '
        )
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def process_rds(self):
        """Mapping result with matched scRNA.
        """
        Mapping.run_mapping(
            self.rds[0], self.contig_file, self.sample, self.outdir, self.assign_file[0]
            )

        meta = pd.read_csv(f'{self.outdir}/{self.sample}_meta.csv')
        metaTB = meta[meta['CellTypes'].isin(self.celltype_set)]
        mappedmeta = meta[meta['Class']=='T/BCR']
        mappedmetaTB = mappedmeta[mappedmeta['CellTypes'].isin(self.celltype_set)]
        
        self.add_metric(
            'Total Cell Number in Matched transcriptome',
            meta.shape[0],
            show=False)
        self.add_metric(
            'Cell Number Successfully Mapped to transcriptome',
            mappedmeta.shape[0],
            show=False)
        self.add_metric(
            'T/B cell Number in Matched transcriptome',
            metaTB.shape[0],
            show=False)
        self.add_metric(
            'Cell Number Successfully Mapped to T/B cell in transcriptome',
            mappedmetaTB.shape[0],
            show=False)

    @utils.add_log
    def process_h5ad(self):
        contig_file = pd.read_csv(self.contig_file)
        assembled_cells = set(contig_file.barcode)
        print(self.h5ad)
        h5ad = sc.read_h5ad(self.h5ad[0])
        h5ad.obs['status'] = 'None'
        h5ad.obs.loc[h5ad.obs.index.isin(assembled_cells), "status"] = "T/BCR"
        
        sc.pl.tsne(h5ad, color='status', add_outline=True, legend_loc='on data',
                legend_fontsize=12, legend_fontoutline=2,frameon=False,
                title='Assembled Cells', palette=['C7','C3'])
        plt.savefig(f"{self.outdir}/{self.sample}_mapping.pdf", dpi=300)

    def run(self):
        if self.rds and self.assign_file:
            self.process_rds()
        if self.h5ad:
            self.process_h5ad()


def mapping(args):
    step_name = 'mapping'
    mapping_obj = Mapping(args, step_name)
    mapping_obj.run()


def get_opts_mapping(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
        parser.add_argument('--match_out', help='assemble result from match directory', required=True)
    return parser