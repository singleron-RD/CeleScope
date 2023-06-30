# Quick start

## Installation
1. Clone repo
```
git clone https://github.com/singleron-RD/CeleScope.git
```

2. Create conda environment and install conda packages. 
If you have conda installed, it is recommended to use [mamba](https://github.com/mamba-org/mamba) (which is a faster replacement for Conda):
```
conda install mamba
cd CeleScope
mamba create -n celescope -y --file conda_pkgs.txt
```
Or you can use [micromamba](https://mamba.readthedocs.io/en/latest/installation.html)
```
curl micro.mamba.pm/install.sh | bash
cd CeleScope
micromamba create -n celescope -y --file conda_pkgs.txt
```


3. Install celescope

Make sure you have activated the `celescope` conda environment before running `pip install celescope`. 
```
conda activate celescope
pip install celescope
```

## Command-line interface

CeleScope contains interfaces `multi_{assay}` to generate pipeline scripts for all assays. Assays can be one of:
|assay|data|kit|
|---|------|--------------|
|rna|single-cell RNA-Seq|GEXSCOPE<sup>TM</sup>|
|dynaseq|single-cell dynamic RNA-Seq|DynaSCOPE<sup>TM</sup>|
|tag|single-cell sample multiplexing|CLindex<sup>TM</sup>|
|vdj|single-cell VDJ|GEXSCOPE<sup>TM</sup> IR|
|flv_CR|single-cell full length VDJ|sCircle<sup>TM</sup>|
|flv_trust4|single-cell full length VDJ|sCircle<sup>TM</sup>|
|capture_virus|single-cell virus|FocuSCOPE<sup>TM</sup> mRNA × EBV|
|snp|single-cell variant|FocuSCOPE<sup>TM</sup>|
|fusion|single-cell fusion|FocuSCOPE<sup>TM</sup>|
|sweetseq|single-cell glycosylation|ProMoSCOPE<sup>TM</sup>|


Run `multi_{assay} -h` to display available arguments. For example：
```
$ multi_rna -h

usage: rna multi-samples [-h] --mapfile MAPFILE [--mod {sjm,shell}] [--queue QUEUE] [--rm_files] [--steps_run STEPS_RUN]
...
```

## Usage Example

Take single-cell rna as an example:

1. Make a rna genomeDir

### Homo sapiens

```
mkdir hs_ensembl_99
cd hs_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz

conda activate celescope
celescope utils mkgtf Homo_sapiens.GRCh38.99.gtf Homo_sapiens.GRCh38.99.filtered.gtf
celescope rna mkref \
 --genome_name Homo_sapiens_ensembl_99_filtered \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.99.filtered.gtf
```

### Mus musculus

```
mkdir mmu_ensembl_99
cd mmu_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz 
gunzip Mus_musculus.GRCm38.99.gtf.gz

conda activate celescope
celescope utils mkgtf Mus_musculus.GRCm38.99.gtf Mus_musculus.GRCm38.99.filtered.gtf

celescope rna mkref \
 --genome_name Mus_musculus_ensembl_99_filtered \
 --fasta Mus_musculus.GRCm38.dna.primary_assembly.fa \
 --gtf Mus_musculus.GRCm38.99.filtered.gtf
```

2. Generate scripts for each sample

Under your working directory, write a shell script `run.sh` as

```
conda activate celescope
multi_rna\
	--mapfile ./rna.mapfile\
	--genomeDir {path to hs_ensembl_99 or mmu_ensembl_99}\
	--thread 8\
	--mod shell
```
`--mapfile` Required.  Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

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
- `sweetseq` Required, matched_dir.
 
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

`--genomeDir` Required. The path of the genome directory after running `celescope rna mkref`.

`--thread` Threads to use. The recommended setting is 8, and the maximum should not exceed 20.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

3. Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## Usage

- [rna](./assay/multi_rna.md)
- [dynaseq](./assay/multi_dynaseq.md)
- [tag](./assay/multi_tag.md)
- [vdj](./assay/multi_vdj.md)
- [flv_CR](./assay/multi_flv_CR.md)
- [flv_trust4](./assay/multi_flv_trust4.md)
- [capture_virus](./assay/multi_capture_virus.md)
- [snp](./assay/multi_snp.md)
- [fusion](./assay/multi_fusion.md)
- [citeseq](assay/multi_citeseq.md)
- [sweetseq](assay/multi_citeseq.md)


## Test scripts and data
https://github.com/singleron-RD/celescope_test_script

 
