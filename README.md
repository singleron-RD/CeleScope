
# CeleScope

## Requirements

- conda
- linux
- minimum 32GB RAM(to run STAR aligner)

## Installation

```
git clone https://github.com/zhouyiqi91/CeleScope.git
cd CeleScope
sh build.sh
```

## Reference genome 

### Homo sapiens

```
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

mkdir -p references/Homo_sapiens/Ensembl/GRCh38
gzip -c -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.fa
gzip -c -d Homo_sapiens.GRCh38.99.gtf.gz > references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.gtf

conda activate celescope

gtfToGenePred -genePredExt -geneNameAsName2 references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.gtf /dev/stdout | \
    awk '{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.refFlat

STAR \
    --runMode genomeGenerate \
    --runThreadN 6 \
    --genomeDir references/Homo_sapiens/Ensembl/GRCh38 \
    --genomeFastaFiles references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.fa \
    --sjdbGTFfile references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.gtf \
    --sjdbOverhang 100
```

### Mus musculus

```
wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

mkdir -p references/Mus_musculus/Ensembl/GRCm38
gzip -c -d Mus_musculus.GRCm38.dna.primary_assembly.fa.gz > references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.fa
gzip -c -d Mus_musculus.GRCm38.99.gtf.gz > references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.99.gtf

conda activate celescope

gtfToGenePred -genePredExt -geneNameAsName2 references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.99.gtf /dev/stdout | \
    awk '{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.99.refFlat

STAR \
    --runMode genomeGenerate \
    --runThreadN 6 \
    --genomeDir references/Mus_musculus/Ensembl/GRCm38 \
    --genomeFastaFiles references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.fa \
    --sjdbGTFfile references/Mus_musculus/Ensembl/GRCm38/Mus_musculus.GRCm38.99.gtf \
    --sjdbOverhang 100
```

## Usage

### Single cell RNA-Seq

```
conda activate celescope
celescope rna run\   
 --fq1 ./data/R2005212_L1_1.fq.gz\
 --fq2 ./data/R2005212_L1_2.fq.gz\
 --chemistry scopeV2.0.1\
 --genomeDir /SGR/references/Homo_sapiens/Ensembl/GRCh38\
 --sample R2005212\
 --thread 4
```

`--fq1` gzipped FASTQ read 1 file path  
`--fq2` gzipped FASTQ read 2 file path  
`--chemistry` chemistry version, choices=['scopeV2.0.0', 'scopeV2.0.1', 'scopeV2.1.0', 'scopeV2.1.1']  
`--genomeDir` reference genome directory path  
`--sample` sample name  
`--thread` number of threads


### Single Cell Multiplexing

```
conda activate celescope
celescope smk run\   
 --fq1 {fq1.gz}\
 --fq2 {fq2.gz}\
 --chemistry scopeV2.0.1\
 --SMK_barcode_fasta {SMK.fasta}\
 --match_dir {match_dir}\
 --dim 2\
 --combine_cluster {combine_cluster.tsv}
```

`--sample_name` output sample prefix  
`--match_dir` matched scRNA-Seq directory after running CeleScope  
`--SMK_barcode_fasta` SMK tag fasta file  
`--dim` SMK dimension
`--combine_cluster` conbine cluster tsv file; first column:original cluster number; second column:combined cluster number  

sample file

```
$cat SMK.fasta
>SMK1
GTTGTCAAGATGCTACCGTTCAGAGATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG
>SMK2
GTTGTCAAGATGCTACCGTTCAGAGTGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG
>SMK3
GTTGTCAAGATGCTACCGTTCAGAGCGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT
>SMK4
GTTGTCAAGATGCTACCGTTCAGAGATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT

$cat combine_cluster.tsv 
1	1
2	2
3	2
4	2
5	3
```



