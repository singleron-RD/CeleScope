STEPS = [
    'mkref',
    'sample',
    'barcode',
    'cutadapt',
    'star',
    'prep',
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


