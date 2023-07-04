__ASSAY__ = 'capture_rna'

STEPS = [
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
