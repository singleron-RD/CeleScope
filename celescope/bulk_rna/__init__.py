STEPS = [
    'mkref',
    'sample',
    'barcode',
    'cutadapt',
    'star',
    'prep_map',
    "featureCounts",
    "count"]
__ASSAY__ = 'bulk_rna'

REMOVE_FROM_MULTI = {
    'mkref',
    'prep_map',
}

IMPORT_DICT = {
    'mkref': 'celescope.rna',
    'star': 'celescope.rna'
}
