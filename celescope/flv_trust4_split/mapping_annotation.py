import celescope.flv_trust4.mapping_annotation as tools_anno
from celescope.tools import utils


class Mapping_annotation(tools_anno.Mapping_annotation):
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


@utils.add_log
def mapping_annotation(args):
    with Mapping_annotation(args, display_title="V(D)J Annotation") as runner:
        runner.run()


def get_opts_mapping_annotation(parser, sub_program):
    tools_anno.get_opts_mapping_annotation(parser, sub_program)