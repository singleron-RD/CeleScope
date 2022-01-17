import os
from collections import OrderedDict

__VERSION__ = "1.7.1"
__version__ = __VERSION__

ASSAY_LIST = [
    "rna",
    'vdj',
    'tag',
    'dynaseq',
    'snp',
    'capture_virus',
    'fusion',
    'hla',
    'capture_rna',
    'citeseq',
    'tcr_fl',
]

ROOT_PATH = os.path.dirname(__file__)

RELEASED_ASSAYS = ['rna', 'vdj', 'tag', 'dynaseq', 'snp', 'capture_virus']

# argument help
HELP_DICT = {
    'match_dir': 'Match celescope scRNA-Seq directory.',
    'gene_list': 'Required. Gene list file, one gene symbol per line. Only results of these genes are reported.',
    'genomeDir': 'Required. Genome directory after running `celescope rna mkref`.',
    'thread': 'Thread to use.',
    'debug': 'If this argument is used, celescope may output addtional file for debugging.',
    'fasta': 'Required. Genome fasta file. Use absolute path or relative path to `genomeDir`.',
    'outdir': 'Output directory.',
    'matrix_dir': 'Match celescope scRNA-Seq matrix directory.',
    'panel': 'The prefix of bed file in `celescope/data/snp/panel/`, such as `lung_1`.',
    'virus_genomeDir': 'Required. Virus genome directory after running `celescope capture_virus mkref`.',
    'threshold_method': 'One of [otsu, auto, hard, none].'
}

# report metrics help
HELP_INFO_DICT = {
    'matched_barcode_number': {
        'display': 'Number of Matched Cells',
        'info': 'cell barcode number of matched scRNA-Seq sample',
    }
}
