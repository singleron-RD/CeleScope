##Features
- Convert barcodes to 10X format.
Output
- `02.convert/barcode_correspond.txt` Recording barcodes correspondence.(if --method sgr)
- `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads as cellranger input.
- `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads as cellranger input.
- `02.convert/{sample}_S1_L001_R3_001.fastq.gz` New R3 reads as cellranger input
- read1, barcode, read2, sample index are associated with R1, R2, R3, I1 respectively.
    R1: Read 1
    R2: Dual index i5 read(10x Barcode)
    R3: Read 2
## Arguments
`--soft_path` soft path for cellranger.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq1` R2 read file.

`--fq2` R2 read file.

`--fq3` R2 read file.

