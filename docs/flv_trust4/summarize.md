## Features
- CDR3 filtering: contain stop condon, length <=5, etc..

- If barcode A's two chains CDR3s are identical to another barcode B, and A's chain abundance is significantly lower than B's, filter A.

- If `--target_cell_barcode` is provided, the UMI counts of all contigs originated from target cells are multiplied by a weight(default: 6.0) to better distinguish signal from background noise. `--target_cell_barcode` comes from the cell type annotation results of the RNA library.

- Cell-calling is similar to the rna cell-calling algorithm.

## Output
- `04.summarize/clonetypes.tsv` High-level description for each clonotype.

- `04.summarize/{sample}_all_contig.csv` High-level and detailed annotation for each contig.

- `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.

- `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of {sample}_all_contig.csv.

- `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
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

`--fq2` Barcode R2 reads.

`--assemble_out` Result of  assemble dirctory.

`--match_dir` Match scRNA-seq directory.

