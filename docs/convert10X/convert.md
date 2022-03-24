## Features

- Format barcodes and UMIs.

## Output        
- `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

- `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads.

- `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads.
## Arguments
`--split_R2` whether split r2.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq2` R2 read file.

