import pandas as pd
import glob
from celescope.tools import utils
from celescope.tools.step import s_common
from celescope.fl_vdj_CR.VDJ_Mixin import VDJ_Mixin, get_opts_VDJ_Mixin


class Mapping(VDJ_Mixin):
    """
    ## Features

    - Assembled T/B cells Mapping with SC-RNA barcodes.

    ## Output
    - `{sample}_assign.png` Auto-assigned umap plot in scRNA-Seq library.

    - `{sample}_cluster_umap.png` Cluster umap plot in scRNA-Seq library.

    - `{sample}_umapplot.png` Umap plot of assembled and un-assembled barcodes in scRNA-Seq library.

    - `{sample}_distribution.txt` Number of assembled barcodes in every clusters in scRNA-Seq library.

    """
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        
        self.seqtype = args.seqtype
        self.match_dir = args.match_dir

        try:
            self.rds = glob.glob(f'{self.match_dir}/06.analysis/*.rds')[0]
            self.assign_file = glob.glob(f'{self.match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')[0]
        except IndexError:
            pass
  
        self.contig = glob.glob(f'{self.outdir}/../05.match/match_contigs.csv')[0]

        if self.seqtype == 'TCR':
            self.Celltype = {'T_cells','NKT_cells','T cells','NK T cells','Tcells'}
            self._name = "Tcells"
        elif self.seqtype == 'BCR':
            self.Celltype = {'Plasma_cells','B_cells','Mature_B_cell', 'Plasma cells', 'B cells','Bcells'}
            self._name = "Bcells"

    @utils.add_log
    def process(self):

        self.run_mapping()
        meta = pd.read_csv(glob.glob(f'{self.outdir}/{self.sample}_meta.csv')[0])
        metaTB = meta[meta['CellTypes'].isin(self.Celltype)]
        mappedmeta = meta[meta['Class']=='T/BCR']
        mappedmetaTB = mappedmeta[mappedmeta['CellTypes'].isin(self.Celltype)]
        
        Transcriptome_cell_number = meta.shape[0]
        TB_cell_number = metaTB.shape[0]
        Mapped_Transcriptome_cell_number = mappedmeta.shape[0]
        Mapped_TB_cell_number = mappedmetaTB.shape[0]

        self.add_metric(
            name='Cell Number in Matched transcriptome',
            value=Transcriptome_cell_number,
            help_info="the number of barcodes considered as cell-associated in matched transcriptome."
            )
        self.add_metric(
            name='Cell Number Successfully Mapped to transcriptome',
            value=Mapped_Transcriptome_cell_number,
            total=Transcriptome_cell_number,
            help_info="Cells with at least one productive chain successfully mapped to transcriptome."
            )
        self.add_metric(
            name=f'{self._name} Number in Matched transcriptome',
            value=TB_cell_number,
            help_info=f"the number of barcodes auto-assigned as {self._name} in matched transcriptome."
            )
        self.add_metric(
            name=f'Cell Number Successfully Mapped to {self._name} in transcriptome',
            value=Mapped_TB_cell_number,
            total=TB_cell_number,
            help_info=f"{self._name} with at least one productive chain successfully mapped to transcriptome."
            )

    def run(self):
        try:
            if self.rds and self.assign_file:
                self.process()
                self._write_stat()
        except AttributeError:
            pass
            # print("rds file and type file do not exist" + "\n" )

def mapping(args):
    step_name = 'mapping'
    mapping_obj = Mapping(args, step_name)
    mapping_obj.run()


def get_opts_mapping(parser, sub_program):
    get_opts_VDJ_Mixin(parser)
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
    return parser