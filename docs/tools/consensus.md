## Features
- Consensus all the reads of the same (barcode, UMI) combinations into one read(UMI).

## Output
- `{sample}_consensus.fq` Consensus fastq.


## Arguments
`--threshold` Default 0.5. Valid base threshold.

`--not_consensus` Skip the consensus step.

`--min_consensus_read` Minimum number of reads to support a base.

`--fq` Required. Fastq file.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

