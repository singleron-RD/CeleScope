__STEPS__ = [
    'mkref',
    'sample',
    'barcode',
    'cutadapt',
    'consensus',
    "star_virus",
    "count_capture_virus",
    "analysis_capture_virus",
    "featureCounts_capture_virus",
    "count_capture_virus_mtx"
]
__ASSAY__ = 'capture_virus'

IMPORT_DICT = {
    'star_virus': 'celescope.rna_virus',
}
