__STEPS__ = [
    'sample',
    'barcode',
    'cutadapt',
    'STAR',
    "featureCounts",
    "count",
    'analysis']
__ASSAY__ = 'rna'

# m: memory 
# x: thread
RESOURCE = {
    'sample': {'m':1, 'x':1},
    'barcode': {'m':5, 'x':1},
    'cutadapt': {'m':5, 'x':1},
    'STAR': {'m':30, 'x':1},
}
