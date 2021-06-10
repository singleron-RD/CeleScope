import os

__VERSION__ = "1.3.1"
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
