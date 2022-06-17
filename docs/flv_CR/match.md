## Features

- Assembled results match with sc-RNA library.

- Generate matched VDJ-annotation metrics, clonetypes table and bar-plot of clonetypes distribution in html.

## Output

- `matched_contig_annotations.csv` High-level annotations of each high-confidence contigs from matched cell-associated barcodes.

- `matched_contig.fasta` High-confidence contig sequences annotated in the matched_contig_annotations.csv.

- `matched_productive_contig_annotations.csv` Annotations of each productive contigs from matched cell-associated barcodes. This is a subset of matched_contig_annotations.csv.

- `matched_productive_contig.fasta` Productive contig sequences annotated in the matched_productive_contig_annotations.csv.

- `clonotypes.csv` High-level descriptions of each clonotype where barcodes match with scRNA-Seq.
## Arguments
`--seqtype` TCR or BCR.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` scRNA-seq match directory.

`--summarize_out` assemble result in SGR barcode from summarize directory.

