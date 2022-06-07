import glob

import celescope.flv_trust4_split.mapping_annotation as split_anno
from celescope.flv_trust4.summarize import Summarize


class Mapping_annotation(split_anno.Mapping_annotation):
    """
    ## Features

    - Assembled T/B cells mapping to transcriptome.
    - Generate VDJ annotation info, clonetypes table and bar-plot of clonetypes distribution in html.

    ## Output
    - `05.mapping_annotation/{sample}_assign.png` Umap plot of Auto-assigned celltype in transcriptome.

    - `05.mapping_annotation/{sample}_cluster_umap.png` Umap plot of Cluster in transcriptome.

    - `05.mapping_annotation/{sample}_umapplot.png` Umap plot of assembled barcodes marked as read color.

    - `05.mapping_annotation/{sample}_distribution.txt` Number of assembled barcodes in every clusters.

    """
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.seqtype = args.seqtype
        self.match_dir = args.match_dir
        self.chains, self.paired_groups = Summarize._parse_seqtype(self.seqtype)
        self.contig_file = args.contig_file

        try:
            self.rds = glob.glob(f'{self.match_dir}/06.analysis/*.rds')[0]
            self.assign_file = glob.glob(f'{self.match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')[0]
        except IndexError:
            pass
  
        self.contig = glob.glob(f'{self.outdir}/../03.summarize/{self.sample}_filtered_contig.csv')[0]

        if self.seqtype == 'TCR':
            self.Celltype = {'T_cells','NKT_cells','T cells','NK T cells','Tcells'}
            self._name = "Tcells"
        elif self.seqtype == 'BCR':
            self.Celltype = {'Plasma_cells','B_cells','Mature_B_cell','Plasma cells','B cells','Bcells'}
            self._name = "Bcells"


def mapping_annotation(args):
    with Mapping_annotation(args, display_title="V(D)J Annotation") as runner:
        runner.run()


def get_opts_mapping_annotation(parser, sub_program):
    split_anno.get_opts_mapping_annotation(parser, sub_program)