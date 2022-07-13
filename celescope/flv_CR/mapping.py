import pandas as pd
import glob
import subprocess
import os
import celescope
from celescope.tools import utils
from celescope.tools.step import s_common, Step


TOOLS_DIR = os.path.dirname(celescope.tools.__file__)

CELL_TYPE_DICT = {
    'TCR':['T_cells', 'NKT_cells', 'T cells', 'NK T cells', 'Tcells'],
    'BCR':['Plasma_cells', 'B_cells', 'Mature_B_cell', 'Plasma cells', 'B cells', 'Bcells'],
    }

class Mapping(Step):
    """
    ## Features

    - Output assembled T/B cells mapping to transcriptome if rds and auto-assign info exist in match directory.

    ## Output
    - `05.annotation/{sample}_assign.png` Umap plot of Auto-assigned celltype in transcriptome.

    - `05.annotation/{sample}_cluster_umap.png` Umap plot of Cluster in transcriptome.

    - `05.annotation/{sample}_umapplot.png` Umap plot of assembled barcodes marked as read color.
    
    - `05.annotation/{sample}_distribution.txt` Number of assembled barcodes in every clusters.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        
        self.seqtype = args.seqtype
        self.match_dir = args.match_dir
        if self.match_dir != 'None':
            self.contig_file = glob.glob(f'{args.match_out}/matched_contig_annotations.csv')[0]
        self.celltype_set = CELL_TYPE_DICT[self.seqtype]

        try:
            self.rds = glob.glob(f'{self.match_dir}/06.analysis/*.rds')[0]
            self.assign_file = glob.glob(f'{self.match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')[0]
        except IndexError:
            pass

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
    def process(self):
        """Mapping result with matched scRNA.
        """
        Mapping.run_mapping(
            self.rds, self.contig_file, self.sample, self.outdir, self.assign_file
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

    def run(self):
        try:
            if self.rds and self.assign_file:
                self.process()
        except AttributeError:
            print("rds file and type file do not exist" + "\n" )
        except ZeroDivisionError:
            print("Not found auto-assigned T/B cell in matched sc-RNA")

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