
# CeleScope
GEXSCOPE Single Cell Analysis Pipelines  
[中文文档](https://github.com/zhouyiqi91/CeleScope/wiki)

## Requirements

- conda
- git
- minimum 32GB RAM(to run STAR aligner)

## Installation


1. Clone repo
`git clone https://github.com/singleron-RD/CeleScope.git`

2. Install conda packages
```
cd CeleScope
conda create -n celescope
conda activate celescope
conda install --file conda_pkgs.txt --channel conda-forge --channel bioconda --channel r --channel imperial-college-research-computing
```

3. Install celescope
```
pip install celescope
# if you are in china, you can use pypi mirror to accelerate downloading
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple celescope
```

4. Install Beta version(optional)
```
# if you want to use Beta version of celescope
python setup.py install
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
 --genomeDir /SGR/references/Homo_sapiens/Ensembl/GRCh38\
 --sample R2005212\
 --chemistry auto\
 --thread 4\
```

`--fq1` Required. Gzipped FASTQ read 1 file path.

`--fq2` Required. Gzipped FASTQ read 2 file path.

`--genomeDir` Required. Reference genome directory.  

`--sample` Required. Sample name. 

`--chemistry` Chemistry version, default=auto. 

`--thread` Number of threads to use, default=1.

### Single Cell VDJ

```
conda activate celescope
celescope vdj run\   
 --fq1 {vdj fq1.gz}\
 --fq2 {vdj fq2.gz}\
 --type {TCR or BCR}\
 --sample {sample name}\
 --chemistry auto\
 --thread {thread}\
 --match_dir {match_dir}\
```  

`--fq1` Required. Gzipped FASTQ read 1 file path.

`--fq2` Required. Gzipped FASTQ read 2 file path.

`--type` Required. TCR or BCR.   

`--sample` Required. Sample name. 

`--chemistry` Chemistry version, default=auto. 

`--thread` Number of threads to use, default=1.

`--match_dir` Optional. Matched scRNA-Seq directory after running CeleScope.  


