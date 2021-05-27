__STEPS__ = [
    'mkref',
    'sample', 'barcode', 'cutadapt', 'consensus', 'star', 'featureCounts', 
    'target_metrics', 'snpCalling', 'analysis_snp'
]
__ASSAY__ = 'snp'
IMPORT_DICT = {
    'star': 'celescope.rna'
}
