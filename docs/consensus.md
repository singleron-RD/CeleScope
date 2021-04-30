# consensus

## Features
- Consensus all the reads of the same (barcode, UMI) combinations into one read(UMI).

## Input
- R2 Fastq.

## Output
- Consensus fastq.

## Parameters

`--fq` Required. R2 reads from step cutadapt.

`--threshold` Default `0.5`. Valid base threshold.

`--not_consensus` Skip the consensus step.

## Metrics

- UMI Counts: Total UMIs number after consensus.

- Mean UMI Length: Mean UMI length after consensus.

- Ambiguous Base Counts: Total ambiguous bases number.
