STEPS = [
    "mkref",
    "sample",
    "barcode",
    "cutadapt",
    "consensus",
    "mapping_vdj",
    "count_vdj",
]
__ASSAY__ = "bulk_vdj"

IMPORT_DICT = {
    "mkref": "celescope.vdj",
}
