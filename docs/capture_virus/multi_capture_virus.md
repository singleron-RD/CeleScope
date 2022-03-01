## Usage
1. Make a virus genomeDir using `celescope capture_virus mkref`

2. Generate shell scripts
```
multi_capture_virus \
--mapfile {mapfile} \
--virus_genomeDir {virus_genomeDir} \
--not_consensus \
--allowNoPolyT \
--thread 4 \
--mod shell
```

3. Run shell scripts
```
sh ./shell/{sample}.sh
```
## Features
### mkref
- Create a virus genome reference directory.


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

### consensus
- Consensus all the reads of the same (barcode, UMI) combinations into one read(UMI). It will go through the sequence residue by residue and 
count up the number of each type of residue (ie. A or G or T or C for DNA) in all sequences in the
alignment. If the following conditions are met, the consensus sequence will be the most common residue in the alignment:
1. the percentage of the most common residue type > threshold(default: 0.5);
2. most common residue reads >= min_consensus_read;
otherwise an ambiguous character(N) will be added.


### star_virus
- Map reads to the viral genome using STAR.


### count_virus

- Taking the bam file as input, count the number of UMIs and reads mapped to the viral genome.


### filter_virus
- Filter background virus UMIs base on a UMI threshold.
There are three methods to determine the UMI threshold:
    - 'auto' : Using a method similar to cell calling method.
    - 'otsu' : UMI counts are first log 2 transformed and then the threshold is determined by [Otsu's method](https://en.wikipedia.org/wiki/Otsu%27s_method)
    - 'hard' : Using User provided UMI threshold.


## Output files
### mkref

- STAR genome index files
- Genome config file


### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### cutadapt
- `cutadapt.log` Cutadapt output log file.
- `{sample}_clean_2.fq.gz` R2 reads file without adapters.

### consensus
- `{sample}_consensus.fq` Fastq file after consensus.

### star_virus
- `{sample}_virus_Aligned.sortedByCoord.out.bam` : Aligned BAM sorted by coordinate.

### count_virus
- {sample_raw_read_count.json} : barcode - UMI - raw_reads_count

### filter_virus
- {sample}_filtered_UMI_tsne.csv : Filtered virus UMI counts of each cell barcode.

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

`--threshold` Default 0.5. Valid base threshold.

`--not_consensus` Skip the consensus step.

`--min_consensus_read` Minimum number of reads to support a base.

`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--outFilterMatchNmin` Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.

`--out_unmapped` Output unmapped reads.

`--STAR_param` Other STAR parameters.

`--outFilterMultimapNmax` Default `1`. How many places are allowed to match a read at most.

`--starMem` Default `30`. Maximum memory that STAR can use.

`--virus_genomeDir` Virus genome dir.

`--min_query_length` Minimum query length.

`--not_correct_UMI` Perform UMI correction.

`--read_threshold_method` method to find read threshold. UMIs with `support reads` < `read threshold` are filtered.

`--read_hard_threshold` int, use together with `--read_threshold_method hard`.

`--umi_threshold_method` method to find UMI threshold. Cell barcode with `UMI` < `UMI threshold` are considered negative.

`--umi_hard_threshold` int, use together with `--umi_threshold_method hard`.

