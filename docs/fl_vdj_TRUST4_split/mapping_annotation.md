## Features

- TCR/BCR full length assembly results.

## Output
- `04.summarize/clonetypes.tsv` High-level descriptions of each clonotype.
- `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
- `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
- `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig_annotations.csv.
- `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
- `04.summarize/{sample}_one_chain_contig.csv`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.csv.
- `04.summarize/{sample}_one_chain_contig.fasta`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.fasta.
## Arguments
`--seqtype` TCR or BCR.

`--coef` coef for auto filter.

`--diffuseFrac` If cell A's two chains CDR3s are identical to another cell B, and A's chain abundance is significantly lower than B's, filter A.

`--expected_target_cell_num` Expected T or B cell number. If `--target_cell_barcode` is provided, this argument is ignored.

`--target_cell_barcode` Barcode of target cells. It is a plain text file with one barcode per line.

`--target_weight` UMIs of the target cells are multiplied by this factor. Only used when `--target_cell_barcode` is provided.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` scRNA-seq match directory.

`--trust_report` Filtered trust report,Filter Nonfunctional CDR3 and CDR3 sequences containing N.

`--barcode_report` Filtered barcode report of trust4 which is related to diffuseFrac option.

`--contig_file` original contig annotation file.

