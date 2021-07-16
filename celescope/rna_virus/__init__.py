__STEPS__ = [
    'sample',
    'barcode',
    'cutadapt',
    'star',
    "star_virus",
    "featureCounts",
    "count",
    "count_virus",
    'analysis_rna_virus',
]
__ASSAY__ = 'rna_virus'
IMPORT_DICT = {
    'star': 'celescope.rna'
}
