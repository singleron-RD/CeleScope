STEPS = [
    "sample",
    "starsolo",
    "cells",
    "split",
    "conversion",
    "substitution",
    "replacement",
]

__ASSAY__ = "bulk_dynaseq"

REMOVE_FROM_MULTI = {
    "cells",
}

IMPORT_DICT = {
    "starsolo": "celescope.bulk_rna",
    "split": "celescope.bulk_rna",
}
