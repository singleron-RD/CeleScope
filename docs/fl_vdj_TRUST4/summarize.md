## Features

- TCR/BCR full length assembly results.

## Output
- `04.summarize/clonetypes.csv` High-level descriptions of each clonotype.
- `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
- `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
- `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig.csv.
- `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
## Arguments
`--seqtype` TCR or BCR.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--full_len_assembly` full len assembly file.

`--assign_out` read assignment results.

`--filter_report` Filtered trust report.

`--barcode_filter_report` Filtered barcode report.

`--cutadapted_fq` Cutadapt R2 reads.

`--assembled_fa` Read used for assembly.

