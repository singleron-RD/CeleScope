# Introduction

V(D)J repertoires Analysis

# Usage example

- Install celescope [celescope](https://github.com/singleron-RD/CeleScope/blob/master/docs/installation.md)

- Download and unpack cellranger soft and reference file.

```
Soft
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-vdj/cellranger-6.1.2.tar.gz?Expires=1646072261&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC12ZGovY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDYwNzIyNjF9fX1dfQ__&Signature=Z-2m906CV5Rb1snIAga-QDSXYSZ8cNqCj1EECGP4uloU3qH~uCMH42MHf4TNnDL2zAsKA7cXsCsQYz0A9yJdNh7dfRT8ohpuAzASFx5Pj-bkqfw4p2tql55IIaPN0zqxyUuyZ9sfKl5qTQX82LoVolRpiBUL8dF9nr~bA2P1gJZ~xg1QssS7icR5MmTzvKKS5NYkezG8vWaTiEdXU0nuKI2ciZSX5GOMeIRW-YYR7mJwHmBbTVxe0o-uBuUtqor0Y98jdIv8Z~dwMjujRjrEShdCGNixTSonGzeS2~9CXqWquCJIOolqFFkcFHgXkD7ZWNfSXWbTxuF57rCsub98pA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

Reference: human and mouse
wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0.tar.gz

tar -xzvf cellranger-6.1.2.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0.tar.gz
```
- Generate scripts for each sample
Under your working directory, write a shell script run.sh as
```
conda activate celescope
multi_vdj_full_len \
    --mapfile  ./test.mapfile \
    --chemistry flv \
    --mem 10 \
    --thread 8 \
    --allowNoLinker \
    --seqtype TCR \
    --ref_path "/SGRNJ/Database/script/soft/cellranger/vdj_ref/6.0.0/hs/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0" \
    --soft_path "/SGRNJ03/randd/cjj/soft/cellranger/cellranger-6.1.2/cellranger" \
```
```
--chemistry flv
--mem Memory(G) for assembly. Default 10.
--thread number of multi-threaded for assembly.
--seqtype Choose from ['TCR', 'BCR']. Required.
--ref_path Absolute path of reference.
--soft_path Absolute path of soft.
```
After you sh run.sh, a shell directory containing {sample}.sh files will be generated.
- Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the ./shell/{sample}.sh must be run under the working directory(You shouldn't run them under the shell directory)
