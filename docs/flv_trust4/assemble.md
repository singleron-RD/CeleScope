## Features

- TRUST4 does not use multi-processing when assembling. By default, the candidate reads are split into 4 chunks to speed up.

- Keep only full-length contigs.

## Output
- `03.assemble/assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.

- `03.assemble/assemble/{sample}_assign.out` Read assignment results.

- `03.assemble/assemble/{sample}_assembled_reads.fa` Assembled raw reads sequence.

- `03.assemble/assemble/{sample}_annotate.fa` Assembled annotated contig sequences info.

- `03.assemble/assemble/{sample}_report.tsv` Record assembled CDR3 types, read count and proportion.

- `03.assemble/assemble/{sample}_filter_report.tsv` Filter non-functional CDR3.

- `03.assemble/assemble/{sample}_barcode_report.tsv` Record chain information for each barcode.

- `03.assemble/assemble/{sample}_barcode_filter_report.tsv` If --diffuseFrac provided. Filter low abundance barcode.
## Arguments
`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--candidate_fq` Candidate fastq file from mapping step.

`--not_split` do not split reads into chunks.

`--ref` reference name.

`--seqtype` TCR/BCR seq data.

`--barcodeRange` Barcode range in fq1, INT INT CHAR.

`--umiRange` UMI range in fq1, INT INT CHAR.

