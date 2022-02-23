## Features

- Summarize contig and clonetypes infomation.

## Output
- `filtered_contig_annotations.csv' High-level annotations of each high-confidence, cellular contig.

- `filtered_contig.fasta` High-confidence contig sequences in cell barcodes.

- `clonotypes.csv` High-level descriptions of each clonotype.

- `match_contigs.csv` Consider barcodes match scRNA-Seq library in filtered_contig_annotations.csv.

- `match_contig.fasta` Consider barcodes match scRNA-Seq library in filtered_contig.fasta.

- `match_clonotypes.csv` Consider barcodes match scRNA-Seq library in clonotypes.csv.


## Arguments
`--not_split_R2` whether split r2

`--seqtype` TCR or BCR

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--barcode_dic` 10X barcode correspond sgr barcode

