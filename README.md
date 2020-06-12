[![issues](https://img.shields.io/github/issues/zhouyiqi91/scope_tools_1.0)](https://github.com/zhouyiqi91/scope_tools_1.0/issues)

# Scope-tools

## Requirements

- conda
- linux
- minimum 32GB RAM(to run STAR aligner)

## Installation

```
git clone https://github.com/SingleronBio/scope_tools_1.0.git
conda env create -f scope1.0.yml
```

## Reference genome 

### Homo sapiens

```
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

mkdir -p references/Homo_sapiens/Ensembl/GRCh38
gzip -c -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.fa
gzip -c -d Homo_sapiens.GRCh38.99.gtf.gz > references/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.99.gtf

conda activate scope1.0

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

conda activate scope1.0

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

### Example

```
conda activate scope1.0
python3 /SGRNJ/Database/script/pipe/scope_tools_1.0/tools/scope.py run  
 --fq1 ./data/R2005212_L1_1.fq.gz\
 --fq2 ./data/R2005212_L1_2.fq.gz\
 --bcType scope\
 --genomeDir /SGRNJ/Public/Database/genome/homo_mus\
 --refFlat /SGRNJ/Public/Database/genome/homo_mus/homo_mus.refFlat\
 --annot /SGRNJ/Public/Database/genome/homo_mus/homo_mus.gtf\
 --outdir R2005212\
 --sample R2005212
```

`--fq1` gzipped FASTQ read 1 file path  
`--fq2` gzipped FASTQ read 2 file path  
`--bcType` barcode and UMI structure type,defalut=scope  
`--genomeDir` reference genome directory path  
`--refFlat` refFlat file path  
`--annot` annotation file path  
`--outdir` output directory path  
`--sample` sample ID  


