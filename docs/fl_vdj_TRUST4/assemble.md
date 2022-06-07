## Features

- Assemble TCR/BCR seq data.

## Output
- `03.assemble/{sample}_R1.fq` record barcode and umi info for assembly
- `03.assemble/{sample}_R1.fq` record sequence and quality info for assembly
- `03.assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.
- `03.assemble/{sample}_assign.out` read assignment results.
- `03.assemble/{sample}_assembled_reads.fa` Assembled raw reads.
- `03.assemble/{sample}_annot.fa` Assembled annotated contig sequences.
- `03.assemble/{sample}_full_len.fa` Assembled full length contig sequences.
- `03.assemble/{sample}_report.tsv` Record assembled CDR3 types, read count and proportion of read count.
- `03.assemble/{sample}_filter_report.tsv` Filter nonfunctional CDR3
- `03.assemble/{sample}_barcode_report.tsv` Record chain information for each barcode.
- `03.assemble/{sample}_barcode_filter_report.tsv` Filter low abundance cell.
## Arguments
`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--cutadapted_fq` cutadapted_fq.

`--match_dir` Match scRNA-seq directory.

`--species` Species name and version.

`--seqtype` TCR/BCR seq data.

`--match_previous_assemble` whether match reads with sc-RNA before assemble.

