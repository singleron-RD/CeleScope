## Download and unpack cellranger soft and reference file.
```
wget -O cellranger-atac-2.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-2.1.0.tar.gz?Expires=1659105655&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hdGFjL2NlbGxyYW5nZXItYXRhYy0yLjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTkxMDU2NTV9fX1dfQ__&Signature=ZzSn3fB9aI4GR9S4uaUOSaH3aKV5aJG2JgPTIQIV5Uka3GYhjc8QSTU~Eb4osKlLo8pghC4ze0PqpwOwUW6UGbhaX~eSoj-Vvei~kSxX3f0psyJ6Kbi4yscWe48RFIRaZkL91pFdApjEskhvDVh2a5eKTpBkt3dolSLncmyoc~7JT-D3dDymYvkZiM2KcfPHBi2lkBp5kxsYVmEzwH3btE0EczIjPCJwZbXHn7CU11XxyrTkp0gUU4Yp2QzglzggJ1kPveSlaUxeZgv~YVs9d-VQ36yVhgu~Yum~tcDZ~1hyr1aR9uigllni8g90uTdKfEaiC6idJ3kh2k9mlWh6iQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

Reference: human, mouse, human-and-mouse
wget https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0.tar.gz

tar -xzvf cellranger-atac-2.1.0.tar.gz
tar -xzvf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
tar -xzvf refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz
tar -xzvf refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0.tar.gz
```

## Usage

```
conda activate celescope
multi_flv_CR \
    --mapfile ./test.mapfile \
    --customized \
    --ref_path "/SGRNJ06/randd/USER/cjj/cr-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0" \
    --soft_path "/SGRNJ06/randd/USER/cjj/cr-atac/cellranger-atac-2.1.0/" \
    --mod shell
```
## Features
### convert

- Convert barcodes to 10X format.

Output

- `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

- `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads as cellranger input.

- `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads as cellranger input.

- `02.convert/{sample}_S1_L001_R3_001.fastq.gz` New R3 reads as cellranger input

- read1, barcode, read2, sample index are associated with R1, R2, R3, I1 respectively.
    R1: Read 1
    R2: Dual index i5 read(10x Barcode)
    R3: Read 2

### atac

- Run cellranger-atac

## Output files
## Arguments
`--mapfile` Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  
4th column: The 4th column has different meaning for each assay. The single cell rna directory after running CeleScope is called `matched_dir`.

- `rna` Optional, forced cell number.
- `vdj` Required, matched_dir.
- `tag` Required, matched_dir.
- `dynaseq` Optional, forced cell number.
- `snp` Required, matched_dir.
- `capture_virus` Required, matched_dir.
- `fusion` Required, matched_dir.
- `citeseq` Required, matched_dir.
- `flv_CR` Required, matched_dir.
- `flv_trust4` Required, matched_dir.
 
5th column:
- `dynaseq` Required, background snp file.

Example

Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.
```
$cat ./my.mapfile
fastq_prefix1	fastq_dir1	sample1
fastq_prefix2	fastq_dir2	sample1
fastq_prefix3	fastq_dir1	sample2

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
```

`--mod` Which type of script to generate, `sjm` or `shell`.

`--queue` Only works if the `--mod` selects `sjm`.

`--rm_files` Remove redundant fastq and bam files after running.

`--steps_run` Steps to run. Multiple Steps are separated by comma. For example, if you only want to run `barcode` and `cutadapt`, 
use `--steps_run barcode,cutadapt`.

`--outdir` Output directory.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--chemistry` chemistry version.

`--gzip` Output gzipped fastq files.

`--bulk_seq` bulk_seq or not.

`--ref_path` reference path for cellranger-atac.

`--soft_path` soft path for cellranger-atac.

`--other_param` Other cellranger-atac count parameters.

`--mem` memory(G).

