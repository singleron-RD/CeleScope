## Features

- Assembled T/B cells mapping to transcriptome.
- Generate VDJ annotation info, clonetypes table and bar-plot of clonetypes distribution in html.

## Output
- `05.mapping_annotation/{sample}_assign.png` Umap plot of Auto-assigned celltype in transcriptome.

- `05.mapping_annotation/{sample}_cluster_umap.png` Umap plot of Cluster in transcriptome.

- `05.mapping_annotation/{sample}_umapplot.png` Umap plot of assembled barcodes marked as read color.

- `05.mapping_annotation/{sample}_distribution.txt` Number of assembled barcodes in every clusters.
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

