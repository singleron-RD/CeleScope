# Quick start

CeleScope contains interfaces `multi_{assay}` to generate pipeline scripts for all assays. Assays can be one of:

- rna
- vdj
- tag
- dynaseq
- snp

Run `multi_{assay} -h` for help.


## Usage Example

Take Single-cell rna as an example:

1. Generate scripts for each sample

Under your working directory, write a shell script `run.sh` as

```
conda activate celescope
multi_rna\
	--mapfile ./rna.mapfile\
	--genomeDir /SGRNJ/Public/Database/genome/homo_mus\
	--thread 8\
	--mod shell
```
`--mapfile` Required. See below on how to write a mapfile.

`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--thread` The recommended setting is 8, and the maximum should not exceed 20.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

2. You can start your analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

3. See [multi_rna.md](./rna/multi_rna.md) for all available arguments.

## Uasge

- [multi_rna.md](./rna/multi_rna.md)
- [multi_vdj.md](./vdj/multi_vdj.md)
- [multi_tag.md](./tag/multi_tag.md)
- [multi_dynaseq.md](./dynaseq/multi_dynaseq.md)
- [multi_snp.md](./snp/multi_snp.md)

## How to write mapfile

Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  
4th column: The 4th column has different meaning for each assay. The single cell rna directory after running CeleScope is called `matched_dir`.
- `rna` Optional, forced cell number.
- `vdj` Optional, matched_dir.
- `tag` Required, matched_dir.
- `dynaseq` Optional, forced cell number.
- `snp` Required, matched_dir.

### Example

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

## Test data
https://github.com/singleron-RD/celescope_test_script


 
