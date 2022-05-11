## Features
- Align R2 reads to IGMT(http://www.imgt.org/) database sequences with mixcr.

## Output
- `{sample}_consensus.fasta` Fasta file after UMI consensus.

- `{sample}_UMI_count_unfiltered.tsv` UMI reading for each (barcode, chain, VJ_pair) combination.

- `{sample}_UMI_count_filtered.tsv` For each (barcode, chain) combination, only the record with the 
most VJ_pair UMI reads is kept.

- `{sample}_align.txt` Result report.

- `{sample}_alignments.txt` The alignment result of each UMI/read.
## Arguments
`--type` TCR or BCR.

`--species` Default `hs`. `hs`(human) or `mmu`(mouse).

`--not_consensus` Input fastq is not consensused.

`--fq` Required. Input fastq file.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

