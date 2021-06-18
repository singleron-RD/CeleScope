import os

__VERSION__ = "1.3.2b1"
__version__ = __VERSION__

ASSAY_DICT = {
    "rna": "Single-cell rna",
    "rna_virus": "Single Cell RNA-Seq Virus",
    'capture_virus': "Single Cell Capture Virus",
    'fusion': "Single Cell Fusion Gene",
    'vdj': "Single-cell vdj",
    'mut': 'Single Cell Insertion and Deletion',
    'hla': 'Single Cell HLA',
    'capture_rna': 'Single Cell Capture RNA',
    'snp': 'Single-cell variant',
    'tag': 'Single-cell tag',
    'citeseq': 'Single Cell CITE-Seq',
    'tcr_fl': 'Single Cell full length TCR',
}

ROOT_PATH = os.path.dirname(__file__)

RELEASED_ASSAYS = ['rna', 'vdj', 'tag', ]

HELP_DICT = {
    'match_dir': 'Match celescope scRNA-Seq directory.',
    'gene_list': 'Gene list file, one gene symbol per line. Only results of these genes are reported.',
    'genomeDir': 'Genome directory after running `mkref`.',
}
