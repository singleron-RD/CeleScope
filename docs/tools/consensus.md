## Features
- Consensus all the reads of the same (barcode, UMI) combinations into one read(UMI).

## Output
- `{sample}_consensus.fq` Consensus fastq.


## Arguments
`--threshold` Default 0.5. Valid base threshold.

`--not_consensus` Skip the consensus step.

`--fq` Required. Fastq file.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

