##Features

- Summarize contig and clonetypes infomation.

- Generate Mapped reads and Cells metrics, barcode rank plot, clonotypes table and clonotypes frequency barplot in html.

Output
- `filtered_contig_annotations.csv` High-level annotations of each high-confidence, cellular contig.

- `filtered_contig.fasta` High-confidence contig sequences in cell barcodes.

- `clonotypes.csv` High-level descriptions of each clonotype.

- `match_contigs.csv` Consider barcodes match scRNA-Seq library in filtered_contig_annotations.csv.

- `match_contig.fasta` Consider barcodes match scRNA-Seq library in filtered_contig.fasta.

- `match_clonotypes.csv` Consider barcodes match scRNA-Seq library in clonotypes.csv.
## Arguments
`--seqtype` TCR or BCR.

`--not_split_R2` not split R2 reads.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--barcode_dict` 10X barcode correspond sgr barcode.

`--assemble_out` assemble result.

`--annotation_out` annotation result.

`--match_out` match result.

