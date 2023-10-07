## Usage
1. Make a vdj reference

### Homo sapiens

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

### Mus musculus

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

2. Generate scripts for each sample

Under your working directory, write a shell script `run.sh` as

```
multi_bulk_vdj \
    --mapfile ./vdj.mapfile \
    --ref_path {path to Homo sapiens or Mus musculus} \
    --species {human or mouse} \
    --type TCR \
    --min_consensus_read 2 \
    --thread 8 \
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

`--ref_path` Required. The path of the reference directory after running `celescope vdj mkref`.

`--species` Required. Human or Mouse.

`--type` Required. TCR or BCR.

`--thread` Threads to use. The recommended setting is 8, and the maximum should not exceed 20.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

After you `sh run.sh`, a `shell` directory containing `{sample}.sh` files will be generated.

3. Start the analysis by running:
```
sh ./shell/{sample}.sh
```
Note that the `./shell/{sample}.sh` must be run under the working directory(You shouldn't run them under the `shell` directory)

## Main output
- `count_vdj/{sample}_clonetypes.csv` The frequency and proportion of each cdr3 amino acid sequence.
- `count_vdj/{sample}_clonetypes_by_nt.csv` The frequency and proportion of each cdr3 nucleotide sequence.
- `count_vdj/{sample}_mapping_metrics.tsv` Mapping metrics of each index.
- `count_vdj/{sample}_filtered_annotations.csv` Annotation for each chain.