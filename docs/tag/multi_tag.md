`multi_tag` is used to analyze Single-cell Multiplexing data generated with CLindex<sup>TM</sup> Sample Multiplexing kits. 
Before running `multi_tag`, you need to run scRNA-Seq data with CeleScope first.

## Usage

```
multi_tag \
    --mapfile ./tag.mapfile\
    --barcode_fasta ./tag_barcode.fasta\
    --fq_pattern L25C15\
    --split_matrix\
    --mod shell
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

### mapping_tag
- Align R2 reads to the tag barcode fasta.


### count_tag
- Assign tag to each cell barcode and summarize.


### analysis_tag
- Combine scRNA-Seq clustering infromation with tag assignment.

### split_tag
- Split scRNA-Seq fastq according to tag assignment.


## Output files
### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### cutadapt
- `cutadapt.log` Cutadapt output log file.
- `{sample}_clean_2.fq.gz` R2 reads file without adapters.

### mapping_tag

- `{sample}_read_count.tsv` tab-delimited text file with 4 columns.

    `barcode` cell barcode  
    `tag_name`  tag name in barcode_fasta  
    `UMI`   UMI sequence  
    `read_count` read count per UMI  

### count_tag

- `{sample}_umi_tag.tsv` 

    `first column` cell barcode  
    `last column`  assigned tag  
    `columns between first and last` UMI count for each tag 

- `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation

- `{sample}_cluster_count.tsv` cell barcode number assigned to *undeterminded*, *multiplet* and *each tag*

### split_tag
- `matrix/` Matrix files of each tag.(Optional)
- `fastq/` Fastq files of each tag.(Optional)

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

`--fq_pattern` Required. R2 read pattern. The number after the letter represents the number of bases.         
`L` linker(common sequences)  
`C` tag barcode.

`--barcode_fasta` Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < len(tag barcode) / 10 + 1. 
If no such tag exists, the read is classified as invalid.

You can find the barcode fasta file under `celescope/data/Clindex`
```
>CLindex_TAG_1
CGTGTTAGGGCCGAT
>CLindex_TAG_2
GAGTGGTTGCGCCAT
>CLindex_TAG_3
AAGTTGCCAAGGGCC
>CLindex_TAG_4
TAAGAGCCCGGCAAG
>CLindex_TAG_5
TGACCTGCTTCACGC
>CLindex_TAG_6
GAGACCCGTGGAATC
>CLindex_TAG_7
GTTATGCGACCGCGA
>CLindex_TAG_8
ATACGCAGGGTCCGA
>CLindex_TAG_9
AGCGGCATTTGGGAC
>CLindex_TAG_10
TCGCCAGCCAAGTCT
>CLindex_TAG_11
ACCAATGGCGCATGG
>CLindex_TAG_12
TCCTCCTAGCAACCC
>CLindex_TAG_13
GGCCGATACTTCAGC
>CLindex_TAG_14
CCGTTCGACTTGGTG
>CLindex_TAG_15
CGCAAGACACTCCAC
>CLindex_TAG_16
CTGCAACAAGGTCGC
```

`--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.

`--UMI_min` Default='auto'. Minimum UMI threshold. Cell barcodes with valid UMI < UMI_min are classified as *undeterminded*.

`--dim` Default=1. Tag dimentions. Usually we use 1-dimentional tag.

`--SNR_min` Default='auto'. Minimum signal-to-noise ratio. 
Cell barcodes with UMI >=UMI_min and SNR < SNR_min are classified as *multiplet*.

`--combine_cluster` Conbine cluster tsv file.

`--coefficient` Default=0.1. If `SNR_min` is 'auto', minimum signal-to-noise ratio is calulated as 
`SNR_min = max(median(SNRs) * coefficient, 2)`. 
Smaller `coefficient` will cause less *multiplet* in the tag assignment.

`--split_fastq` If used, will split scRNA-Seq fastq file according to tag assignment.

`--split_matrix` If used, will split scRNA-Seq matrix file according to tag assignment.

`--split_vdj` If used, will split scRNA-Seq vdj count file according to tag assignment.

`--vdj_dir` Match celescope vdj directory. Required when --split_vdj is specified.

