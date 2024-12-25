# barcode
PATTERN_DICT = {
    "auto": None,
    "scopeV1": "C12U8T18",
    "scopeV2.0.0": "C8L16C8L16C8U8T18",
    "scopeV2.0.1": "C8L16C8L16C8L1U8T18",
    "scopeV2.1.0": "C8L16C8L16C8U12T18",
    "scopeV2.1.1": "C8L16C8L16C8L1U12T18",
    "scopeV2.2.1": "C8L16C8L16C8L1U12T18",
    "scopeV3.0.1": "C9L16C9L16C9L1U12T18",
    # flv_rna is actually U9L16, but the last 10 bases can not be sequenced accurately.
    "flv_rna": "C8L16C8L16C8U9L6",
    "flv": "U9C8L16C8L16C8",
    "bulk_vdj": "L18C6U16",
    "bulk_rna": "C9U12",
    "customized": None,
    "rna_5p": "U8C9L4C9L4C9",
    "rna_3p": "C9L4C9L4C9U8",
    "5p3p-1": "C9C9C9U8",  # converted
    # temp
    "rna_5p.1": "C8L4U9C8L4C8",
    "rna_3p.1": "C8L4C8U9L4C8",
    "5p3p-2": "C8C8C8U9",  # converted
    "rna_5p.2": "C9L6U9C9L4C9",
    "rna_3p.2": "C9L4C9U9L6C9",
    "5p3p-3": "C9C9C9U9",  # converted
}


# count
OUTS_DIR = "outs"
RAW_MATRIX_DIR_SUFFIX = "raw"
FILTERED_MATRIX_DIR_SUFFIX = "filtered"
MATRIX_FILE_NAME = "matrix.mtx.gz"
FEATURE_FILE_NAME = "features.tsv.gz"
BARCODE_FILE_NAME = "barcodes.tsv.gz"
STAR_BAM_SUFFIX = "Aligned.out.bam"
TAG_BAM_SUFFIX = "aligned_posSorted_addTag.bam"
STARSOLO_BAM_SUFFIX = "Aligned.sortedByCoord.out.bam"
COUNTS_FILE_NAME = "counts.tsv"

# mkref
GENOME_CONFIG = "celescope_genome.config"
