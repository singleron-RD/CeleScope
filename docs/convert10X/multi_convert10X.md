## Installation

1. Clone repo
```
git clone -b convert_10X https://github.com/singleron-RD/CeleScope.git
```

2. Create conda environment and install conda packages
```
cd CeleScope
conda create -n convert_10X -y --file conda_pkgs.txt
```

Alternatively, you can use [mamba](https://github.com/mamba-org/mamba) to improve speed.
```
conda install mamba
mamba create -n convert_10X -y --file conda_pkgs.txt
```

3. Install celescope

Make sure you have activated the conda environment before running `pip install Celescope`. 
```
conda activate convert_10X
pip install .
```

## Usage
For rna data, create run.sh file as shown below:
```
conda activate convert_10X
multi_convert10X \
    --mapfile  test.mapfile \
    --chemistry flv_rna \
```
For vdj data, create run.sh file as shown below:
```
conda activate convert_10X
multi_convert10X \
    --mapfile  test.mapfile \
    --chemistry flv \
```
Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.
1st column: Fastq file prefix.
2nd column: Fastq file directory path.
3rd column: Sample name, which is the prefix of all output files.
```
rna     /SGRNJ03/randd/cjj/celedev/TESTDATA/testcele/celescope_test_data/rna/fastqs/    test1
```
## Features
### barcode

- Demultiplex barcodes.
- Filter invalid R1 reads, which includes:
    - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
    - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
    - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
    - Low quality reads: low sequencing quality in barcode and UMI regions.


### convert

- Format barcodes and UMIs.


## Output files
### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### convert
- `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

- `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads in 10X format.

- `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads in 10X format.

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

`--chemistry` Predefined (pattern, barcode whitelist, linker whitelist) combinations. `--chemistry auto` can auto-detect scopeV2 mRNA, scopeV3 mRNA, full length VDJ mRNA(flv_rna) and full length VDJ(flv). You need to explicitly use `--chemistry scopeV1` for legacy chemistry scopeV1. `--chemistry customized` is used for user defined combinations that you need to provide `--pattern`, `--whitelist` and `--linker` at the same time.

`--pattern` The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  
- `C`: cell barcode  
- `L`: linker(common sequences)  
- `U`: UMI    
- `T`: poly T.

`--whitelist` Cell barcode whitelist file path, one cell barcode per line.

`--linker` Linker whitelist file path, one linker per line.

`--lowQual` Default 0. Bases in cell barcode and UMI whose phred value are lower than lowQual will be regarded as low-quality bases.

`--lowNum` The maximum allowed lowQual bases in cell barcode and UMI.

`--nopolyT` Outputs R1 reads without polyT.

`--noLinker` Outputs R1 reads without correct linker.

`--allowNoPolyT` Allow valid reads without polyT.

`--allowNoLinker` Allow valid reads without correct linker.

`--gzip` Output gzipped fastq files.

`--output_R1` Output valid R1 reads.

`--split_R2` whether split r2.

