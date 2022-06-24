## Features

- V(D)J results match SC-RNA infomation.

## Output
- `match_contigs.csv` Consider barcodes match scRNA-Seq library in filtered_contig_annotations.csv.

- `match_contig.fasta` Consider barcodes match scRNA-Seq library in filtered_contig.fasta.

- `match_clonotypes.csv` Consider barcodes match scRNA-Seq library in clonotypes.csv.
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

`--match_dir` scRNA-seq match directory.

