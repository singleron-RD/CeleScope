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


### sample
- Generate sample info.
- Add read_ID in fastq file. The format of the read name is `{readId}_{index}`.


### cutadapt
- Trim adapters in R2 reads with cutadapt. Default adapters includes:
    - polyT=A{18}, 18 A bases. 
    - p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Illumina p5 adapter.

### mapping_vdj
- Align R2 reads to IGMT(http://www.imgt.org/) database sequences with blast.


### count_vdj
- Summarize clonetypes infomation.


## Output files
### mkref

- VDJ IMGT reference with index files.

### sample
- `01.barcode/{sample}_2.fq(.gz)`

### cutadapt
- `cutadapt.log` Cutadapt output log file.
- `{sample}_clean_2.fa` R2 reads file without adapters in fasta format.

### mapping_vdj
- `{sample}_airr.tsv` The alignment result of each read.
A tab-delimited file compliant with the AIRR Rearrangement schema(https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

- `{sample}_produtive.tsv` Including all productive chains for each readID.

### count_vdj
- `{sample}_filtered_annotations.csv` Annotations for each CDR3.

- `{sample}_clonetypes.csv` The count and percentage of each CDR3 of reads.

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

`--debug` If this argument is used, celescope may output addtional file for debugging.

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

`--ref_path` reference path for igblast.

`--species` Default human. human or mouse.

`--type` Required. `TCR` or `BCR`.

