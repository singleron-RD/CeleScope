## Features

- V(D)J results match SC-RNA infomation.

## Output
- `match_contigs.csv` Consider barcodes match scRNA-Seq library in filtered_contig_annotations.csv.

- `match_contig.fasta` Consider barcodes match scRNA-Seq library in filtered_contig.fasta.

- `match_clonotypes.csv` Consider barcodes match scRNA-Seq library in clonotypes.csv.


## Arguments
`--seqtype` TCR or BCR

`--remove_3rd_chain` remove IGK or IGL according to umis when a cell has 3 chains at the same time.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--barcode_dic` 10X barcode correspond sgr barcode

`--match_dir` scRNA-seq match directory

