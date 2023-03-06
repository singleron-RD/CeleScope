##Features

- Convert barcodes and UMI to 10X format.

Output

- `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

- `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads as cellranger input.

- `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads as cellranger input.
## Arguments
`--soft_path` soft path for cellranger.

`--tenX_chemistry` 10X chemistry version, V2 or V3 for scRNA, V2 for VDJ.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq2` R2 read file.

