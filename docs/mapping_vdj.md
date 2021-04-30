# mapping_vdj

## Features
- Align R2 reads to IGMT(http://www.imgt.org/) database sequences with mixcr.

## Input
- R2 reads from step cutadapt.

## Output
- `{sample}_consensus.fasta` Fasta file after UMI consensus.

- `{sample}_UMI_count_unfiltered.tsv` UMI reading for each (barcode, chain, VJ_pair) combination.

- `{sample}_UMI_count_filtered1.tsv` For each (barcode, chain) combination, only the record with the most VJ_pair UMI reads is kept.

- `{sample}_align.txt` Result report.

- `{sample}_alignments.txt` The alignment result of each UMI/read.

## Parameters

`--fq` Required. R2 reads from step cutadapt. 

`--type` Required. `TCR` or `BCR`.

`--species` Default `hs`. `hs` or `mmu`.


## Metrics
If '--not_consensus' argument was used, Reads were used instead of UMIs.

- UMIs(or Reads) Mapped to Any VDJ Gene : UMIs(or Reads) Mapped to any germline VDJ gene segments.

- UMIs(or Reads) with CDR3 : UMIs(or Reads) with CDR3 sequence.

- UMIs(or Reads) with Correct CDR3 : some UMIs(or Reads) with CDR3 might have stop codon and these UMIs(or Reads) are classified as incorrect.

- UMIs(or Reads) Mapped Confidently to VJ Gene : UMIs(or Reads) mapped to VJ gene pairs and with correct CDR3.

- UMIs(or Reads) Mapped to IGH : UMIs(or Reads) mapped confidently to IGH chain.

- UMIs(or Reads) Mapped to IGL : UMIs(or Reads) mapped confidently to IGL chain.

- UMIs(or Reads) Mapped to IGK : UMIs(or Reads) mapped confidently to IGK chain.