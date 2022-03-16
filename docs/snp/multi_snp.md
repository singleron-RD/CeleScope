## Usage

### Make a snp reference genomeDir

1. Run `celescope rna mkref`. If you already have a rna genomeDir, you can use it and skip this step.
2. Run `celescope snp mkref` under the rna genomeDir. Check [mkref.md](./mkref.md) for help.

### Install ANNOVAR, download the annotation database and write a annovar config file.
https://annovar.openbioinformatics.org/en/latest/user-guide/download/

```
perl /Public/Software/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar cosmic70 humandb/
```

annovar_config file
```
[ANNOVAR]
dir = /Public/Software/annovar/  
db = /SGRNJ/Database/script/database/annovar/humandb  
buildver = hg38  
protocol = refGene,cosmic70  
operation = g,f  
```

### Run multi_snp
There are two ways to run `multi_snp`

1. Do not perform consensus before alignment and report read count(recommended for data generated with FocuSCOPE kit).

```
multi_snp\
    --mapfile ./test1.mapfile\
    --genomeDir {genomeDir after running celescope snp mkref}\
    --thread 4\
    --mod shell\
    --panel lung_1\
    --annovar_config annovar.config\
    --not_consensus
```

2. Do consensus before alignment and report UMI count. 

```
multi_snp\
    --mapfile ./test1.mapfile\
    --genomeDir {genomeDir after running celescope snp mkref}\
    --thread 4\
    --mod shell\
    --panel lung_1\
    --annovar_config annovar.config\
```
## Features
### mkref
- Create dictionary file and fasta index for gatk SplitNCigarReads.
(https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) 
Need to run `celescope rna mkref` first


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


### star
- Align R2 reads to the reference genome with STAR.
- Collect Metrics with Picard.


### featureCounts
- Assigning uniquely mapped reads to genomic features with FeatureCounts.

### target_metrics
- Filter bam file
    - Filter reads that are not cell-associated.
    - Filter reads that are not mapped to target genes. 

- Collect enrichment metrics.


### variant_calling
- Perform variant calling at single cell level.


### filter_snp
- Filter out `ref` and `alt` alleles that do not have enough reads to support.


### analysis_snp
- Annotate variants with [Annovar](https://annovar.openbioinformatics.org/en/latest/).


## Output files
### mkref
- fasta index
- gatk dictionary file


### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### cutadapt
- `cutadapt.log` Cutadapt output log file.
- `{sample}_clean_2.fq.gz` R2 reads file without adapters.

### consensus
- `{sample}_consensus.fq` Fastq file after consensus.

### star
- `{sample}_Aligned.sortedByCoord.out.bam` BAM file contains Uniquely Mapped Reads.

- `{sample}_SJ.out.tab` SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format.

- `{sample}_Log.out` Main log with a lot of detailed information about the run. 
This is most useful for troubleshooting and debugging.

- `{sample}_Log.progress.out` Report job progress statistics, such as the number of processed reads, 
% of mapped reads etc. It is updated in 1 minute intervals.

- `{sample}_Log.Log.final.out` Summary mapping statistics after mapping job is complete, 
very useful for quality control. The statistics are calculated for each read (single- or paired-end) and 
then summed or averaged over all reads. Note that STAR counts a paired-end read as one read, 
(unlike the samtools agstat/idxstats, which count each mate separately). 
Most of the information is collected about the UNIQUE mappers 
(unlike samtools agstat/idxstats which does not separate unique or multi-mappers). 
Each splicing is counted in the numbers of splices, which would correspond to 
summing the counts in SJ.out.tab. The mismatch/indel error rates are calculated on a per base basis, 
i.e. as total number of mismatches/indels in all unique mappers divided by the total number of mapped bases.

- `{sample}_region.log` Picard CollectRnaSeqMetrics results.

### featureCounts
- `{sample}` Numbers of reads assigned to features (or meta-features).
- `{sample}_summary` Stat info for the overall summrization results, including number of 
successfully assigned reads and number of reads that failed to be assigned due to 
various reasons (these reasons are included in the stat info).
- `{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam` featureCounts output BAM, 
sorted by coordinates；BAM file contains tags as following(Software Version>=1.1.8):
    - CB cell barcode
    - UB UMI
    - GN gene name
    - GX gene id
- `{sample}_name_sorted.bam` featureCounts output BAM, sorted by read name.

### target_metrics
- `filtered.bam` BAM file after filtering. Reads that are not cell-associated or not mapped to target genes are filtered.

### variant_calling
- `{sample}_raw.vcf` Variants are called with bcftools default settings.
- `{sample}_norm.vcf` Indels are left-aligned and normalized. See https://samtools.github.io/bcftools/bcftools.html#norm for more details.

### filter_snp
- `{sample}_test1_filtered.vcf` VCF file after filtering. Alleles read counts that do not have enough reads to support are set to zero. 
Genotypes are changed accordingly.

### analysis_snp
- `{sample}_gt.csv` Genotypes of variants of each cell. Rows are variants and columns are cells.
- `{sample}_variant_ncell.csv` Number of cells with each genotype.
- `{sample}_variant_table.csv` `{sample}_variant_ncell.csv` annotated with COSMIC(https://cancer.sanger.ac.uk/cosmic).

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

`--outFilterMatchNmin` Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.

`--out_unmapped` Output unmapped reads.

`--STAR_param` Other STAR parameters.

`--outFilterMultimapNmax` Default `1`. How many places are allowed to match a read at most.

`--starMem` Default `30`. Maximum memory that STAR can use.

`--gtf_type` Specify feature type in GTF annotation.

`--featureCounts_param` Other featureCounts parameters.

`--gene_list` Required. Gene list file, one gene symbol per line. Only results of these genes are reported. Conflict with `--panel`.

`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--panel` The prefix of bed file in `celescope/data/snp/panel/`, such as `lung_1`. Conflict with `--gene_list`.

`--threshold_method` One of [otsu, auto, hard, none].

`--hard_threshold` int, use together with `--threshold_method hard`.

`--annovar_config` ANNOVAR config file.

