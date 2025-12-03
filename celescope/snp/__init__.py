STEPS = [
    "mkref",
    "sample",
    "starsolo",
    "target_metrics",
    "variant_calling",
    "filter_snp",
    "analysis_snp",
]
__ASSAY__ = "snp"

PANEL = {"lung_1", "blood_1", "CHIP"}

REMOVE_FROM_MULTI = {
    "mkref",
}
