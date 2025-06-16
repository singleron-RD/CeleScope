## Usage

```
multi_flv_trust4 \
    --mapfile ./test.mapfile \
    --species human or mouse \
    --thread 8 \
    --seqtype BCR \
    --mod shell
```


## Arguments
`--mapfile` Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  
4th column: The single cell rna directory after running CeleScope is called `matched_dir`.

Example

Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.
```
$cat ./my.mapfile
fastq_prefix1	fastq_dir1	sample1 sample1_matched_rna
fastq_prefix2	fastq_dir2	sample1 sample1_matched_rna
fastq_prefix3	fastq_dir1	sample2 sample2_matched_rna

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
```

`--species` If the species is not human or mouse, there is no need to provide the `--species` parameter. Instead, provide the `--ref` parameter. The `--ref` can be constructed using the [BuildImgtAnnot.pl](https://github.com/liulab-dfci/TRUST4?tab=readme-ov-file#practical-notes) script provided by TRUST4.  
For example:  
`BuildImgtAnnot.pl Oryctolagus cuniculus > rabbit_IMGT+C.fa`

`--seqtype` TCR or BCR.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 



## Features
### barcode

- Demultiplex barcodes.
- Filter invalid R1 reads, which includes:
    - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
    - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
    - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
    - Low quality reads: low sequencing quality in barcode and UMI regions.


### mapping

- Extract candidate reads to assemble.


### assemble

- TRUST4 does not use multi-processing when assembling. By default, the candidate reads are split into 4 chunks to speed up.

- Keep only full-length contigs.


### summarize
- CDR3 filtering: contain stop condon, length <=5, etc..

- If barcode A's two chains CDR3s are identical to another barcode B, and A's chain abundance is significantly lower than B's, filter A.

- If `--target_cell_barcode` is provided, the UMI counts of all contigs originated from target cells are multiplied by a weight(default: 6.0) to better distinguish signal from background noise. `--target_cell_barcode` comes from the cell type annotation results of the RNA library.

- Cell-calling is similar to the rna cell-calling algorithm.


### annotation

- Output assembled T/B cells mapping to transcriptome if rds and auto-assign info exist in match directory.
- Generate VDJ annotation info, clonetypes table and bar-plot of clonetypes distribution in html.


## Output files
### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### mapping
- `02.mapping/{sample}_bcrtcr.fq` All candidate reads(mapped to any V(D)J genes) sequence.
- `02.mapping/{sample}_bcrtcr_bc.fa` All candidate reads(mapped to any V(D)J genes) barcode.
- `02.mapping/{sample}_bcrtcr_umi.fa` All candidate reads(mapped to any V(D)J genes) umi.
- `02.mapping/{sample}_{chain}.fq` Candidate reads(mapped to {chain}) sequence. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.
- `02.mapping/{sample}_{chain}_bc.fa` Candidate reads(mapped to {chain}) barcode. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.
- `02.mapping/{sample}_{chain}_umi.fa` Candidate reads(mapped to {chain}) umi. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.

### assemble
- `03.assemble/assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.

- `03.assemble/assemble/{sample}_assign.out` Read assignment results.

- `03.assemble/assemble/{sample}_assembled_reads.fa` Assembled raw reads sequence.

- `03.assemble/assemble/{sample}_annotate.fa` Assembled annotated contig sequences info.

- `03.assemble/assemble/{sample}_report.tsv` Record assembled CDR3 types, read count and proportion.

- `03.assemble/assemble/{sample}_filter_report.tsv` Filter non-functional CDR3.

- `03.assemble/assemble/{sample}_barcode_report.tsv` Record chain information for each barcode.

- `03.assemble/assemble/{sample}_barcode_filter_report.tsv` If --diffuseFrac provided. Filter low abundance barcode.

### summarize
- `04.summarize/clonetypes.tsv` High-level description for each clonotype.

- `04.summarize/{sample}_all_contig.csv` High-level and detailed annotation for each contig.

- `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.

- `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of {sample}_all_contig.csv.

- `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.

### annotation
- `05.annotation/{sample}_assign.png` Umap plot of Auto-assigned celltype in transcriptome.
- `05.annotation/{sample}_cluster_umap.png` Umap plot of Cluster in transcriptome.
- `05.annotation/{sample}_umapplot.png` Umap plot of assembled barcodes marked as read color.
- `05.annotation/{sample}_distribution.txt` Number of assembled barcodes in every clusters.
