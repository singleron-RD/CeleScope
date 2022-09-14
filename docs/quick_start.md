# Quick start

## Installation
Document can be found [here](installation.md)

## Command-line interface

CeleScope contains interfaces `multi_{assay}` to generate pipeline scripts for all assays. Assays can be one of:

- rna
- vdj
- tag
- dynaseq
- snp
- capture_virus
- fusion
- citeseq
- flv_CR
- flv_trust4

Run `multi_{assay} -h` for help.

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
	--genomeDir /SGRNJ/Public/Database/genome/homo_mus\
	--thread 8\
	--mod shell
```
`--mapfile` Required. Check [multi_rna.md](./rna/multi_rna.md) for how to write the mapfile.

`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--thread` Threads to use. The recommended setting is 8, and the maximum should not exceed 20.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

3. Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## Usage

- [multi_rna.md](./rna/multi_rna.md)
- [multi_vdj.md](./vdj/multi_vdj.md)
- [multi_tag.md](./tag/multi_tag.md)
- [multi_dynaseq.md](./dynaseq/multi_dynaseq.md)
- [multi_snp.md](./snp/multi_snp.md)
- [multi_capture_virus.md](./capture_virus/multi_capture_virus.md)
- [multi_fusion.md](./fusion/multi_fusion.md)
- [multi_citeseq.md](citeseq/multi_citeseq.md)
- [multi_flv_CR.md](./flv_CR/multi_flv_CR.md)
- [multi_flv_trust4.md](./flv_trust4/multi_flv_trust4.md)

## Test scripts and data
https://github.com/singleron-RD/celescope_test_script

 
