__STEPS__ = [    
    'sample',
    'barcode',
    'cutadapt',
    "STAR_virus",
    "count_capture_virus",
    "analysis_capture_virus",
]
__ASSAY__ = 'capture_virus'

IMPORT_DICT = {
    'STAR_virus': 'celescope.rna_virus',
}