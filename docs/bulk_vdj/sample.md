## Features
- Generate sample info.
- Add read_ID in fastq file. The format of the read name is `{readId}_{index}`.

## Output
- `01.barcode/{sample}_2.fq(.gz)`
## Arguments
`--gzip` Output gzipped fastq files.

`--output_R1` Output valid R1 reads.

`--fq1` R1 fastq file. Multiple files are separated by comma.

`--fq2` R2 fastq file. Multiple files are separated by comma.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

