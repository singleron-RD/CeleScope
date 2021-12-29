__STEPS__ = [
    'sample',
    'barcode',
    'cutadapt',
    'star',
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

# m: memory
# x: thread
RESOURCE = {
    'sample': {'m': 1, 'x': 1},
    'barcode': {'m': 5, 'x': 1},
    'cutadapt': {'m': 5, 'x': 1},
    'star': {'m': 30, 'x': 1},
}
