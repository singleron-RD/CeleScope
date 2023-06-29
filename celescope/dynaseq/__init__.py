STEPS = [
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


