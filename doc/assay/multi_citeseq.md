## Usage
```
multi_citeseq \
    --mapfile ./test.mapfile \
    --barcode_fasta ./TotalSeqA_mouse_barcode.fasta \
    --fq_pattern L21C15\
    --mod shell
```


## Arguments

`--mapfile` Mapfile is a tab-delimited text file with four columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  
4th column: TThe single cell rna directory after running CeleScope.

Example

Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.
```
$cat ./my.mapfile
fastq_prefix1	fastq_dir1	sample1 sample1_scRNA_celescope_dir
fastq_prefix2	fastq_dir2	sample1 sample1_scRNA_celescope_dir
fastq_prefix3	fastq_dir1	sample2 sample2_scRNA_celescope_dir

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz

$ls sample1_scRNA_celescope_dir
00.sample  01.starsolo  02.analysis  outs  sample1_report.html
```

`--barcode_fasta` A fasta file with read name as protein and sequence as protien barcode.

```
$ head ./TotalSeqA_mouse_barcode.fasta
>Ms.CD4
AACAAGACCCTTGAG
>Ms.CD8a
TACCCGTAATAGCGT
>Ms.CD366
ATTGGCACTCAGATG
>Ms.CD279
GAAAGTCAAAGCACT
>Ms.Ly.6C
AAGTCGTGAGGCATG
```

`--fq_pattern`  R2 read pattern. The number after the letter represents the number of bases.
`L` linker(common sequences)  `C` protein barcode.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

## Main Output

`{sample}/03.count_cite/{sample}_citeseq.mtx.gz` Protein UMI matrix in tsv format. Rows are proteins and columns are cell barcodes.
