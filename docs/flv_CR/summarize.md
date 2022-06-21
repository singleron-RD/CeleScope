## Features

- Convert 10X barcode of assemble result back to SGR barcode.

- Generate Productive contigs sequences and annotation files.

- Generate VDJ-annotation metrics in html.

## Output

- `filtered_contig_annotations.csv` High-level annotations of each high-confidence contigs from cell-associated barcodes.

- `filtered_contig.fasta` High-confidence contig sequences annotated in the filtered_contig_annotations.csv.

- `productive_contig_annotations.csv` Annotations of each productive contigs from cell-associated barcodes. This is a subset of filtered_contig_annotations.csv.

- `productive_contig.fasta` Productive contig sequences annotated in the productive_contig_annotations.csv.

- `clonotypes.csv` High-level descriptions of each clonotype.
## Arguments
`--seqtype` TCR or BCR.

`--soft_path` soft path for cellranger.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--barcode_convert_json` json file.

`--assemble_out` directory of cellranger assemble result.

