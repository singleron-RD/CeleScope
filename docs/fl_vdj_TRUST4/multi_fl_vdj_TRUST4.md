## Installation

1. Clone repo
```
git clone https://github.com/singleron-RD/CeleScope.git
```

2. Create conda environment and install conda packages
```
cd CeleScope
conda create -n celescope -y --file conda_pkgs.txt
```

Alternatively, you can use [mamba](https://github.com/mamba-org/mamba) to improve speed.
```
conda install mamba
mamba create -n celescope -y --file conda_pkgs.txt
```

3. Install celescope

Make sure you have activated the conda environment before running `pip install Celescope`. 
```
conda activate celescope
pip install .
```

## Usage

```
conda activate TRUST_dev
    multi_fl_vdj_TRUST4 \
    --mapfile ./test.mapfile \
    --outdir ./ \
    --chemistry flv \
    --allowNoLinker \
    --species GRCm38 \
    --thread 10 \
    --seqtype BCR \
    --match_previous_assemble
```
## Features
### barcode

- Demultiplex barcodes.
- Filter invalid R1 reads, which includes:
    - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
    - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
    - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
    - Low quality reads: low sequencing quality in barcode and UMI regions.


### cutadapt
- Trim adapters in R2 reads with cutadapt. Default adapters includes:
    - polyT=A{18}, 18 A bases. 
    - p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Illumina p5 adapter.

### assemble

- Assemble TCR/BCR seq data.


### summarize

- TCR/BCR full length assembly results.


## Output files
### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### cutadapt
- `cutadapt.log` Cutadapt output log file.
- `{sample}_clean_2.fq.gz` R2 reads file without adapters.

### assemble
- `03.assemble/{sample}_R1.fq` record barcode and umi info for assembly
- `03.assemble/{sample}_R1.fq` record sequence and quality info for assembly
- `03.assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.
- `03.assemble/{sample}_assign.out` read assignment results.
- `03.assemble/{sample}_assembled_reads.fa` Assembled raw reads.
- `03.assemble/{sample}_annot.fa` Assembled annotated contig sequences.
- `03.assemble/{sample}_full_len.fa` Assembled full length contig sequences.
- `03.assemble/{sample}_report.tsv` Record assembled CDR3 types, read count and proportion of read count.
- `03.assemble/{sample}_filter_report.tsv` Filter nonfunctional CDR3
- `03.assemble/{sample}_barcode_report.tsv` Record chain information for each barcode.
- `03.assemble/{sample}_barcode_filter_report.tsv` Filter low abundance cell.

### summarize
- `04.summarize/clonetypes.csv` High-level descriptions of each clonotype.
- `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
- `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
- `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig.csv.
- `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.

## Arguments
`--mapfile` Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  
4th column: The 4th column has different meaning for each assay. The single cell rna directory after running CeleScope is called `matched_dir`.

- `rna` Optional, forced cell number.
- `vdj` Optional, matched_dir.
- `tag` Required, matched_dir.
- `dynaseq` Optional, forced cell number.
- `snp` Required, matched_dir.
- `capture_virus` Required, matched_dir.
- `fl_vdj_CR` Required, matched_dir.
- `fl_vdj_TRUST4` Required, matched_dir.
- `fl_vdj_TRUST4_split` Required, matched_dir.

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

`--chemistry` Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:  
- `auto` Default value. Used for Singleron GEXSCOPE libraries >= scopeV2 and automatically detects the combinations.  
- `scopeV1` Used for legacy Singleron GEXSCOPE scopeV1 libraries.  
- `customized` Used for user defined combinations. You need to provide `pattern`, `whitelist` and `linker` at the 
same time.

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

`--output_R1` Output valid R1 reads.

`--gzip` Output gzipped fastq files.

`--adapter_fasta` Addtional adapter fasta file.

`--minimum_length` Default `20`. Discard processed reads that are shorter than LENGTH.

`--nextseq_trim` Default `20`. Quality trimming of reads using two-color chemistry (NextSeq). 
Some Illumina instruments use a two-color chemistry to encode the four bases. 
This includes the NextSeq and the NovaSeq. 
In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
However, dark cycles also occur when sequencing “falls off” the end of the fragment.
The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.

`--overlap` Default `10`. Since Cutadapt allows partial matches between the read and the adapter sequence,
short matches can occur by chance, leading to erroneously trimmed bases. 
For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter. 
To reduce the number of falsely trimmed bases, the alignment algorithm requires that 
at least {overlap} bases match between adapter and read.

`--insert` Default `150`. Read2 insert length.

`--cutadapt_param` Other cutadapt parameters. For example, --cutadapt_param "-g AAA".

`--species` Species name and version.

`--match_previous_assemble` whether match reads with sc-RNA before assemble.

`--seqtype` TCR or BCR.

