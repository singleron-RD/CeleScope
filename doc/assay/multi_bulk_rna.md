## Usage
1. Make a rna genomeDir

Please refer to [rna](multi_rna.md) assay.

2. Generate scripts for each sample

Under your working directory, write a shell script `run.sh` as

```
multi_bulk_rna\
	--mapfile ./rna.mapfile\
	--genomeDir {path to hs_ensembl_99 or mmu_ensembl_99}\
	--thread 8\
	--mod shell
```
`--mapfile` Required.  Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  

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

1. Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)




## Features
### sample
- Generate sample info.

### barcode
- Demultiplex barcodes.

### cutadapt
- Trim poly A tails and user-provided adapters in R2 reads with cutadapt.

### star
- Align R2 reads to the reference genome with STAR.

### featureCounts
- Assigning uniquely mapped reads to genomic features with FeatureCounts.

### count
- Generate expression matrix.
- Well statistic.


## Output files

### barcode
- `{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of
    the read name is `{barcode}_{UMI}_{read ID}`.

### cutadapt
- `{sample}_cutadapt.json` Cutadapt output log file.
- `{sample}_clean_2.fq` R2 reads file without adapters in fastq format.

### star
- `{sample}_Aligned.out.bam` BAM file contains Uniquely Mapped Reads.
- `{sample}_Log.final.out` Summary mapping statistics.

### featureCounts
- `{sample}` Numbers of reads assigned to features (or meta-features).
- `{sample}_summary` Stat info for the overall summrization results.
- `{sample}_aligned_posSorted_addTag.bam` featureCounts output BAM, sorted by coordinates, contains tags.

### count
- `{sample}_raw_feature_bc_matrix` The expression matrix of all detected barcodes in [Matrix Market Exchange Formats](
    https://math.nist.gov/MatrixMarket/formats.html).
- `{sample}_counts.txt` 4 columns:
	- Well: barcode sequence
	- UMI: UMI count for each barcode
	- read: read count of each barcode
	- gene: gene count for each barcode
- `{sample}_counts_report.txt` 4 columns as `{sample}_counts.txt`. Filtered by UMI/gene/read cutoff (default UMI>=500).

## Arguments
`--umi_cutoff` If the UMI number exceeds the threshold, it is considered a valid well and reported. (default 500)

`--gene_cutoff`  If the gene number exceeds the threshold, it is considered a valid well and reported. (default 0)

`--read_cutoff`  If the read number exceeds the threshold, it is considered a valid well and reported. (default 0)
