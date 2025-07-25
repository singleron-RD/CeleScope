STEPS = [
    "mkref",
    "sample",
    "barcode",
    "consensus",
    "mapping_vdj",
    "count_vdj",
]
__ASSAY__ = "vdj"

CHAINS = {
    "TCR": ["TRA", "TRB"],
    "BCR": ["IGH", "IGL", "IGK"],
}

PAIRS = {
    "TCR": [["TRA", "TRB"]],
    "BCR": [["IGH", "IGK"], ["IGH", "IGL"]],
}
