import pandas as pd
import glob
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


class Mapping(Step):
    """
    ## Features

    - Output TSNE-plot of Assembled T/B Cells.

    ## Output
    - `06.mapping/{sample}_mapping.pdf` TSNE-plot of Assembled Cells.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        
        self.match_dir = args.match_dir
        if self.match_dir != 'None':
            self.contig_file = glob.glob(f'{args.match_out}/matched_contig_annotations.csv')[0]

        self.h5ad = glob.glob(f'{self.match_dir}/06.analysis/*.h5ad')

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
        if self.h5ad:
            self.process_h5ad()


def mapping(args):
    step_name = 'mapping'
    mapping_obj = Mapping(args, step_name)
    mapping_obj.run()


def get_opts_mapping(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
        parser.add_argument('--match_out', help='assemble result from match directory', required=True)
    return parser