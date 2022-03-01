## Features
- Consensus all the reads of the same (barcode, UMI) combinations into one read(UMI). It will go through the sequence residue by residue and 
count up the number of each type of residue (ie. A or G or T or C for DNA) in all sequences in the
alignment. If the following conditions are met, the consensus sequence will be the most common residue in the alignment:
1. the percentage of the most common residue type > threshold(default: 0.5);
2. most common residue reads >= min_consensus_read;
otherwise an ambiguous character(N) will be added.

## Output
- `{sample}_consensus.fq` Fastq file after consensus.
## Arguments
`--threshold` Default 0.5. Valid base threshold.

`--not_consensus` Skip the consensus step.

`--min_consensus_read` Minimum number of reads to support a base.

`--fq` Required. Fastq file.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

