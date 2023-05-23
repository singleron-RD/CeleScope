## Reference
- Human
```
mkdir -p /genome/vdj/human
cd /genome/vdj/human
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TR{A,B}{V,J}.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBD.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IG{H,K,L}{V,J}.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta
celescope vdj mkref human TR
celescope vdj mkref human IG
```

- Mouse
```
mkdir -p /genome/vdj/mouse
cd /genome/vdj/mouse
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TR{A,B}{V,J}.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TRBD.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IG{H,K,L}{V,J}.fasta
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGHD.fasta
celescope vdj mkref mouse TR
celescope vdj mkref mouse IG
```

## Usage

Please note `multi_vdj` is only used for CDR3 vdj data(not full length). For full length vdj data, please use `multi_flv_CR` or  `multi_flv_trust4`.

```
multi_vdj \
    --mapfile ./vdj.mapfile \
    --ref_path /genome/vdj/human/human_TR \
    --species human \
    --type TCR \
    --thread 8 \
    --mod shell
``` 
## Features
### mkref

- Build Index for IMGT_ref.

- Make sure current directory contains all V,D,J of TRA/TRB or IGH/IGK/IGL reference downloaded from IMGT website.


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


### mapping_vdj
- Align R2 reads to IGMT(http://www.imgt.org/) database sequences with blast.


### count_vdj
- Cell-calling based on barcode-UMI rank.    
- Summarize clonetypes infomation.


## Output files
### mkref

- VDJ IMGT reference with index files.

### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### cutadapt
- `cutadapt.log` Cutadapt output log file.
- `{sample}_clean_2.fq.gz` R2 reads file without adapters.

### consensus
- `{sample}_consensus.fq` Fastq file after consensus.

### mapping_vdj
- `{sample}_airr.tsv` The alignment result of each UMI.
A tab-delimited file compliant with the AIRR Rearrangement schema(https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

- `{sample}_UMI_count_unfiltered.tsv` UMI reading for each (barcode, chain, VJ_pair) combination.

- `{sample}_UMI_count_filtered.tsv` For each (barcode, chain) combination, only the record with the 
most VJ_pair UMI reads is kept.

### count_vdj
- `{sample}_cell_confident.tsv` The clone type of VDJ cell barcode, each chain occupies one line.

- `{sample}_cell_confident_count.tsv` The clone type of VDJ cell barcode, each cell occupies one line.

- `{sample}_clonetypes.tsv` The count and percentage of each clonetypes of VDJ cell barcode.

- `{sample}_match_clonetypes.tsv` When summarize clonetypes, only consider barcodes in the match scRNA-Seq library. 
This file will only be produced when the `match_dir` parameter is provided.

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
- `citeseq` Required, matched_dir.
- `flv_CR` Required, matched_dir.
- `flv_trust4` Required, matched_dir.
- `sweetseq` Required, matched_dir.
 
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

`--use_R3` ATAC libraries use R3 reads instead of R2.

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

`--filterNoPolyT` Filter reads without PolyT.

`--allowNoLinker` Allow valid reads without correct linker.

`--output_R1` Output valid R1 reads.

`--gzip` Output gzipped fastq files.

`--adapter_fasta` Addtional adapter fasta file.

`--minimum_length` Discard processed reads that are shorter than LENGTH.

`--nextseq_trim` Quality trimming of reads using two-color chemistry (NextSeq). 
Some Illumina instruments use a two-color chemistry to encode the four bases. 
This includes the NextSeq and the NovaSeq. 
In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
However, dark cycles also occur when sequencing “falls off” the end of the fragment.
The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.

`--overlap` Since Cutadapt allows partial matches between the read and the adapter sequence,
short matches can occur by chance, leading to erroneously trimmed bases. 
For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter. 
To reduce the number of falsely trimmed bases, the alignment algorithm requires that 
at least {overlap} bases match between adapter and read.

`--insert` Read2 insert length.

`--cutadapt_param` Other cutadapt parameters. For example, --cutadapt_param "-g AAA".

`--threshold` Default 0.5. Valid base threshold.

`--not_consensus` Skip the consensus step.

`--min_consensus_read` Minimum number of reads to support a base.

`--ref_path` reference path for igblast.

`--species` Default human. human or mouse.

`--split_fasta` split fasta file to map to avoid running out of memory.

`--type` Required. `TCR` or `BCR`.

`--UMI_min` minimum number of chain UMI to consider as as cell.

`--BCR_iUMI` Minimum number of UMI of identical receptor type and CDR3 for BCR. 
For each (barcode, chain) combination, only UMI>=iUMI is considered valid.

`--TCR_iUMI` Minimum number of UMI of identical receptor type and CDR3 for TCR. 
For each (barcode, chain) combination, only UMI>=iUMI is considered valid.

`--expected_target_cell_num` Expected T or B cell number. If `--target_cell_barcode` is provided, this argument is ignored.

`--target_cell_barcode` Barcode of target cells. It is a plain text file with one barcode per line. If provided, `--expected_target_cell_num` is ignored.

`--target_weight` UMIs of the target cells are multiplied by this factor. Only used when `--target_cell_barcode` is provided.

