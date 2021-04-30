# consensus

## Features
- Sort fastq by read name.
- Consensus read in name-sorted fastq.

## Input
- Name-sorted fastq.

## Output
- Consensus fastq.

## Parameters

`--fq` Required. R2 reads from step cutadapt.

`--threshold` Default `0.5`. Valid base threshold.

`--not_consensus` Input fastq is not consensus.

## Metrics

- UMI Counts: Total UMIs num.
- Mean UMI Length: Mean UMI length after consensus.
- Ambiguous Base Counts: Total ambiguous bases num.
