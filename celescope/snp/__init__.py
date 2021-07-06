__STEPS__ = [
    'mkref',
    'sample', 'barcode', 'cutadapt', 'consensus', 'star', 'featureCounts',
    'target_metrics', 'variant_calling', 'analysis_snp'
]
__ASSAY__ = 'snp'
IMPORT_DICT = {
    'star': 'celescope.rna'
}
