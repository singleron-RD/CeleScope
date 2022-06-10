## Features

- V(D)J annotation infomation.

## Output
- `filtered_contig_annotations.csv' High-level annotations of each high-confidence, cellular contig.

- `filtered_contig.fasta` High-confidence contig sequences in cell barcodes.

- `clonotypes.csv` High-level descriptions of each clonotype.
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

`--barcode_dict` 10X barcode correspond sgr barcode.

