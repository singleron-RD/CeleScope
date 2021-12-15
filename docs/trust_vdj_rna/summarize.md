## Features

- TCR/BCR full length assembly results.

## Output
- `04.summarize/clonetypes.tsv` High-level descriptions of each clonotype.
- `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
- `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
- `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig_annotations.csv.
- `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
- `04.summarize/{sample}_chain_filtered_contig.csv`Keep the 2 contigs with the highest UMI. This is a subset of filtered_contig.csv.
- `04.summarize/{sample}_chain_filtered_contig.fasta`Keep the 2 contigs with the highest UMI. This is a subset of filtered_contig.fasta.
- `04.summarize/{sample}_one_chain_contig.csv`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.csv.
- `04.summarize/{sample}_one_chain_contig.fasta`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.fasta.
 


## Arguments
`--seqtype` TCR or BCR

`--min_read_count` filter cell by read count number, int type required

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--reads_assignment` File records reads assigned to contigs.

`--fq2` Cutadapt R2 reads.

`--assembled_fa` Read used for assembly

