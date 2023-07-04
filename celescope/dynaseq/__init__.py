STEPS = [
    'sample',
    'barcode',
    'cutadapt',
    'star',
    'prep_map',
    "featureCounts",
    "count",
    'analysis',
    'conversion',
    'substitution',
    'replacement',
    'replace_tsne']

__ASSAY__ = 'dynaseq'

IMPORT_DICT = {
    'star': 'celescope.rna',
    'analysis': 'celescope.rna',
}


REMOVE_FROM_MULTI = {
    'barcode',
    'cutadapt',
    'star',
}
