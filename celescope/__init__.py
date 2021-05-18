import os

__VERSION__ = "1.2.0"
__version__ = __VERSION__

ASSAY_DICT = {
    "rna": "Single Cell RNA-Seq",
    "rna_virus": "Single Cell RNA-Seq Virus",
    'capture_virus': "Single Cell Capture Virus",
    'fusion': "Single Cell Fusion Gene",
    'vdj': "Single Cell VDJ",
    'mut': 'Single Cell Insertion and Deletion',
    'hla': 'Single Cell HLA',
    'capture_rna': 'Single Cell Capture RNA',
    'snp': 'Single Cell Variant',
    'tag': 'Single Cell tag',
    'citeseq': 'Single Cell CITE-Seq',
    'tcr_fl': 'Single Cell full length TCR',
}
ROOT_PATH = os.path.dirname(__file__)
