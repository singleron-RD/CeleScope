from importlib import resources

rna_dict = {
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
        "linker": ["linker1.txt", "linker2.txt"],
    },
    "flv_rna-V2": {
        # 0-3 bases offset at the begining
        "pattern": "L18C9L6C9L6C9U12",
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
        "linker": ["linker1.txt", "linker2.txt", "linker3.txt"],
    },
}

full_length_vdj_dict = {
    "flv": {
        "pattern": "U9C8L16C8L16C8",
        "bc": ["bc.txt", "bc.txt", "bc.txt"],
        "linker": ["linker1.txt", "linker2.txt"],
    },
    "flv-V2": {
        "pattern": "U12C9L6C9L6C9L18",
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
        "linker": ["linker1.txt", "linker2.txt", "linker3.txt"],
    },
}

bulk_vdj_dict = {
    "bulk_vdj": {
        "pattern": "L18C6U16",
        "bc": ["bc.txt"],
    },
}

bulk_rna_dict = {
    "bulk_rna-V1": {
        "pattern": "C9U12",
        "bc": ["bc.txt"],
    },
    "bulk_rna-V2": {
        "pattern": "L9C9U12",
        "bc": ["bc.txt"],
    },
    "bulk_rna-V3": {
        "pattern": "C6U16",
        "bc": ["bc.txt"],
    },
    "bulk_rna-bulk_vdj_match": {
        "pattern": "L18C6U16",
        "bc": ["bc.txt"],
    },
}

space_dict = {
    # 96*96 spots
    "space-ffpe": {
        "pattern": "U12C8L8C8L8",
        "bc": ["bc.txt", "bc.txt"],
        "linker": ["linker1.txt", "linker2.txt"],
    },
    # 144*144 spots
    "space-ffpe-V1.1": {
        "pattern": "U12C8L8C8L8",
        "bc": ["bc.txt", "bc.txt"],
        "linker": ["linker1.txt", "linker2.txt"],
    },
    # fresh frozen
    "space-ff": {
        "pattern": "U12C8L8C8L8L16",
        "bc": ["bc1.txt", "bc2.txt"],
        "linker": ["linker1.txt", "linker2.txt", "linker3.txt"],
    },
}

chemistry_dict = {
    "auto": {
        "pattern": ""  # auto detect
    },
    "customized": {
        "pattern": ""  # user defined
    },
}

for d in [rna_dict, full_length_vdj_dict, bulk_vdj_dict, bulk_rna_dict, space_dict]:
    chemistry_dict.update(d)

chemistry_dir = str(resources.files("celescope.data.chemistry"))
