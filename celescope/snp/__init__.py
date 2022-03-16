__STEPS__ = [
    'mkref',
    'sample', 'barcode', 'cutadapt', 'consensus', 'star', 'featureCounts',
    'target_metrics', 'variant_calling', 'filter_snp', 'analysis_snp'
]
__ASSAY__ = 'snp'
IMPORT_DICT = {
    'star': 'celescope.rna'
}
