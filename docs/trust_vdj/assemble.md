## Features

- Assemble TCR/BCR seq data.

## Output
- `03.assemble/match/{sample}_matched_R1.fq` New R1 reads matched with scRNA-seq
- `03.assemble/match/{sample}_matched_R1.fq` New R2 reads matched with scRNA-seq

- `03.assemble/assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.
- `03.assemble/assemble/{sample}_assign.out` read assignment results.
- `03.assemble/assemble/{sample}_assembled_reads.fa` Assembled raw reads.
- `03.assemble/assemble/{sample}_annot.fa` Assembled annotated contig sequences.
- `03.assemble/assemble/{sample}_full_len.fa` Assembled full length contig sequences.
- `03.assemble/assemble/report.out` Record assembled CDR3 types, read count and proportion of read count.
- `03.assemble/assemble/barcoderep.tsv` Record chain information for each barcode.
- `03.assemble/assemble/barcoderepfl.tsv` Record chain information for each barcode(preliminary filter).
## Arguments
`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq2` R2 reads matched with scRNA-seq.

`--match_dir` Match scRNA-seq directory.

`--species` Species name and version.

`--seqtype` TCR/BCR seq data.

`--barcodeRange` Barcode range in fq1, INT INT CHAR.

`--umiRange` UMI range in fq1, INT INT CHAR.

`--trimLevel` INT: 0: no trim; 1: trim low quality; 2: trim unmatched.

`--UMI_min` UMI number, INT.

