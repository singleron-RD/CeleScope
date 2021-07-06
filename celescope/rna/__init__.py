__STEPS__ = [
    'mkref',
    'sample',
    'barcode',
    'cutadapt',
    'star',
    "featureCounts",
    "count",
    'analysis']
__ASSAY__ = 'rna'

# m: memory
# x: thread
RESOURCE = {
    'sample': {'m': 1, 'x': 1},
    'barcode': {'m': 5, 'x': 1},
    'cutadapt': {'m': 5, 'x': 1},
    'star': {'m': 30, 'x': 1},
}
