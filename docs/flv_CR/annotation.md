## Features

- Convert 10X barcode of assemble result back to SGR barcode.

- Generate VDJ annotation metrics in html.

## Output

- `filtered_contig_annotations.csv' High-level annotations of each high-confidence, cellular contig.

- `filtered_contig.fasta` High-confidence contig sequences in cell barcodes.

- `clonotypes.csv` High-level descriptions of each clonotype.
## Arguments
`--seqtype` TCR or BCR.

`--version` cellranger version.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--barcode_dict` 10X barcode correspond sgr barcode.

`--assemble_out` directory of cellranger assemble result.

