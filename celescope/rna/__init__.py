STEPS = [
    'mkref',
    'sample',
    'barcode',
    'cutadapt',
    'star',
    'prep_map',
    "featureCounts",
    "count",
    'analysis']
__ASSAY__ = 'rna'
REMOVE_FROM_MULTI = {
    'mkref',
    'barcode',
    'cutadapt',
    'star',
}


