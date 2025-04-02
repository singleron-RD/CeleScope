from importlib import resources

chemistry_dict = {
    "auto": {
        "pattern": ""  # auto detect
    },
    "customized": {
        "pattern": ""  # user defined
    },
    "GEXSCOPE-MicroBead": {
        "pattern": "C12U8",
        "bc": [],
    },
    "GEXSCOPE-V1": {
        "pattern": "C8L16C8L16C8L1U12",
        "bc": ["bc.txt", "bc.txt", "bc.txt"],
    },
    "GEXSCOPE-V2": {
        "pattern": "C9L16C9L16C9L1U12",
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
    },
    "GEXSCOPE-V3": {
        # 0-3 bases offset at the begining
        "pattern": "C9L6C9L6C9L1U12",
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
        "linker": ["linker1.txt", "linker2.txt"],
    },
    "flv_rna": {
        "pattern": "C8L16C8L16C8U9L6",
        "bc": ["bc.txt", "bc.txt", "bc.txt"],
    },
    "flv_rna-V2": {
        # 0-3 bases offset at the begining
        "pattern": "L18C9L6C9L6C9U12",
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
        "linker": ["linker1.txt", "linker2.txt", "linker3.txt"],
    },
    "flv": {
        "pattern": "U9C8L16C8L16C8",
        "bc": ["bc.txt", "bc.txt", "bc.txt"],
    },
    "rna_5p-1": {
        "pattern": "U8C9L4C9L4C9",
    },
    "rna_3p-1": {
        "pattern": "C9L4C9L4C9U8",
    },
    "5p3p-1": {
        "pattern": "C9C9C9U8",  # converted
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
    },
    "rna_5p-2": {
        "pattern": "C8L4U9C8L4C8",
    },
    "rna_3p-2": {
        "pattern": "C8L4C8U9L4C8",
    },
    "5p3p-2": {
        "pattern": "C8C8C8U9",  # converted
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
    },
    "rna_5p-3": {
        "pattern": "C9L6U9C9L4C9",
    },
    "rna_3p-3": {
        "pattern": "C9L4C9U9L6C9",
    },
    "5p3p-3": {
        "pattern": "C9C9C9U9",  # converted
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
    },
    "bulk_vdj": {
        "pattern": "L18C6U16",
        "bc": ["bc.txt"],
    },
    "bulk_rna-V1": {
        "pattern": "C9U12",
        "bc": ["bc.txt"],
    },
    "bulk_rna-V2": {
        "pattern": "L9C9U12",
        "bc": ["bc.txt"],
    },
    "bulk_rna-bulk_vdj_match": {
        "pattern": "L18C6U16",
        "bc": ["bc.txt"],
    },
}

chemistry_dir = str(resources.files("celescope.data.chemistry"))
