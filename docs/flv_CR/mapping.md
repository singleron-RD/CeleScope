## Features

- Output assembled T/B cells mapping to transcriptome if rds and auto-assign info exist in match directory.

## Output
- `05.annotation/{sample}_assign.png` Umap plot of Auto-assigned celltype in transcriptome.

- `05.annotation/{sample}_cluster_umap.png` Umap plot of Cluster in transcriptome.

- `05.annotation/{sample}_umapplot.png` Umap plot of assembled barcodes marked as read color.

- `05.annotation/{sample}_distribution.txt` Number of assembled barcodes in every clusters.
## Arguments
`--seqtype` TCR or BCR.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` scRNA-seq match directory.

`--match_out` assemble result from match directory.

