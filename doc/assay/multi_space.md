## Usage
1. Create a genomeDir. 

The genomeDir for spatial data differs from the one used for regular single-cell RNA-seq. Since spatial analysis needs to include some short non-coding RNAs, the GTF file used is unfiltered (the `celescope utils mkgtf` was not applied).

### Homo sapiens

```
mkdir hs_ensembl_110
cd hs_ensembl_110

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz

conda activate celescope
celescope rna mkref \
 --genome_name Homo_sapiens_ensembl_110 \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.110.gtf
```

### Mus musculus

```
mkdir mmu_ensembl_110
cd mmu_ensembl_110

wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.110.gtf.gz

conda activate celescope

celescope rna mkref \
 --genome_name Mus_musculus_ensembl_110 \
 --fasta Mus_musculus.GRCm39.dna.primary_assembly.fa \
 --gtf Mus_musculus.GRCm39.110.gtf 
```

2. Generate scripts for each sample.

Under your working directory, write a shell script `run.sh` as

```
multi_space \
    --mapfile ./sample.mapfile \
    --genomeDir {path to hs_ensembl_110 or mmu_ensembl_110} \
    --thread 16 \
    --mod shell
```
`--mapfile` Required.  Mapfile is a tab-delimited text file with four columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  
4rd column: Spatial directory path.

Example

Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.
```
$cat ./my.mapfile
fastq_prefix1	fastq_dir1	sample1 sample1_spatial_dir
fastq_prefix2	fastq_dir2	sample1 sample1_spatial_dir
fastq_prefix3	fastq_dir1	sample2 sample2_spatial_dir

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
```

**Spatial directory** 
The fourth column of the mapfile is the path to the spatial directory, which can be generated using [AtlasXbrowser](https://github.com/singleron-RD/AtlasXbrowser).

`--genomeDir` Required. The path of the genome directory after running `celescope rna mkref`.

`--thread` It is recommended to use 16 threads. Using more than 20 threads is not advised because  [the mapping speed of STAR saturates at >~20 threads](https://github.com/singleron-RD/CeleScope/issues/197).

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

Start the analysis by running:
```
bash ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## Main output

In the outs folder, there is a `filtered_feature_bc_matrix.h5` file and a `spatial` directory. This format is consistent with 10X Visium output and can be directly loaded into Seurat.

**Seurat v5**

```R
library(Seurat) # v5
obj <- Load10X_Spatial(data.dir = "./outs")
```

**Seurat v4**

Seurat v4 does not recognize the `tissue_positions.parquet` file in the spatial directory and requires the `tissue_positions_list.csv` file. You can first rename the `positions_list.csv` file to `tissue_positions_list.csv`:
```
mv outs/spatial/positions_list.csv outs/spatial/tissue_positions_list.csv
```
then 
```R
library(Seurat) # v4
obj <- Load10X_Spatial(data.dir = "./outs")
```
