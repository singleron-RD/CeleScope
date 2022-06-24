## Features

- Assembled T/B cells Mapping with SC-RNA barcodes.

## Output
- `{sample}_assign.png` Auto-assigned umap plot in scRNA-Seq library.

- `{sample}_cluster_umap.png` Cluster umap plot in scRNA-Seq library.

- `{sample}_umapplot.png` Umap plot of assembled and un-assembled barcodes in scRNA-Seq library.

- `{sample}_distribution.txt` Number of assembled barcodes in every clusters in scRNA-Seq library.
## Arguments
`--species` species.

`--soft` cellranger version.

`--ref_path` reference path for cellranger.

`--soft_path` soft path for cellranger.

`--other_param` Other cellranger parameters.

`--not_split_R2` whether split r2.

`--mem` memory (G).

`--seqtype` TCR or BCR.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` scRNA-seq match directory.

