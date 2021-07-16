__ASSAY__ = 'capture_rna'

__STEPS__ = [
    'sample',
    'barcode',
    'cutadapt',
    'star',
    "featureCounts",
    "count_capture_rna",
    'analysis']

IMPORT_DICT = {
    'star': 'celescope.rna',
    'analysis': 'celescope.rna',
}
