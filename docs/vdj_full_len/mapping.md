## Features

- Assembled T/B cells Mapping with SC-RNA barcodes.

## Output
- `{sample}_assign.png` Auto-assigned umap plot in scRNA-Seq library.

- `{sample}_cluster_umap.png` Cluster umap plot in scRNA-Seq library.

- `{sample}_umapplot.png` Umap plot of assembled and un-assembled barcodes in scRNA-Seq library.

- `{sample}_distribution.txt` Number of assembled barcodes in every clusters in scRNA-Seq library.


## Arguments
`--seqtype` TCR or BCR

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` scRNA-seq match directory

