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
    --thread 10\
    --mod shell\
    --gene_list gene_list.tsv\
    --annovar_config annovar.config\
    --not_consensus
```

2. Do consensus before alignment and report UMI count. 

```
multi_snp\
    --mapfile ./test1.mapfile\
    --genomeDir {genomeDir after running celescope snp mkref}\
    --thread 10\
    --mod shell\
    --gene_list gene_list.tsv\
    --annovar_config annovar.config\
    --min_support_read 1
```


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

`--rm_files` Remove redundant fastq and bam files after running.

`--steps_run` Steps to run. Multiple Steps are separated by comma.

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
- `T`: poly T

`--whitelist` Cell barcode whitelist file path, one cell barcode per line.

`--linker` Linker whitelist file path, one linker per line.

`--lowQual` Default 0. Bases in cell barcode and UMI whose phred value are lower than lowQual will be regarded as low-quality bases.

`--lowNum` The maximum allowed lowQual bases in cell barcode and UMI.

`--nopolyT` Outputs R1 reads without polyT.

`--noLinker` Outputs R1 reads without correct linker.

`--allowNoPolyT` Allow valid reads without polyT.

`--allowNoLinker` Allow valid reads without correct linker.

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

`--threshold` Default 0.5. Valid base threshold.

`--not_consensus` Skip the consensus step.

`--min_consensus_read` Minimum number of reads to support a base.

`--outFilterMatchNmin` Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.

`--out_unmapped` Output unmapped reads.

`--STAR_param` Other STAR parameters.

`--outFilterMultimapNmax` Default `1`. How many places are allowed to match a read at most.

`--starMem` Default `30`. Maximum memory that STAR can use.

`--gtf_type` Specify feature type in GTF annotation

`--featureCounts_param` Other featureCounts parameters

`--gene_list` Required. Gene list file, one gene symbol per line. Only results of these genes are reported.

`--genomeDir` Required. Genome directory after running `mkref`.

`--min_support_read` Minimum number of reads support a variant. If `auto`(default), otsu method will be used to determine this value.

`--annovar_config` ANNOVAR config file.

