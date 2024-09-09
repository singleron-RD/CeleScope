STEPS = [
    "mkref",
    "sample",
    "barcode",
    "cutadapt",
    "consensus",
    "star_virus",
    "count_virus",
    "filter_virus",
    "analysis_virus",
    "featureCounts",
    "count",
]
__ASSAY__ = "capture_virus"

IMPORT_DICT = {
    "star_virus": "celescope.rna_virus",
}
