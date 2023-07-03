## Usage

```
multi_flv_trust4 \
    --mapfile ./test.mapfile \
    --ref hg38 \
    --thread 8 \
    --seqtype BCR \
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

`--gzip` Output gzipped fastq files.

`--output_R1` Output valid R1 reads.

`--not_split` do not split reads into chunks.

`--barcodeRange` Barcode range in fq1, INT INT CHAR.

`--umiRange` UMI range in fq1, INT INT CHAR.

`--ref` reference name.

`--coef` coef for auto filter.

`--diffuseFrac` If cell A's two chains CDR3s are identical to another cell B, and A's chain abundance is significantly lower than B's, filter A.

`--expected_target_cell_num` Expected T or B cell number. If `--target_cell_barcode` is provided, this argument is ignored.

`--target_cell_barcode` Barcode of target cells. Auto or path of plain text file with one barcode per line.

`--target_weight` UMIs of the target cells are multiplied by this factor. Only used when `--target_cell_barcode` is provided.

`--seqtype` TCR or BCR.

