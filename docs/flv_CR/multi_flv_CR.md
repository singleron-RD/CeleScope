## Download and unpack cellranger soft and reference file.
```
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-vdj/cellranger-6.1.2.tar.gz?Expires=1646072261&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC12ZGovY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDYwNzIyNjF9fX1dfQ__&Signature=Z-2m906CV5Rb1snIAga-QDSXYSZ8cNqCj1EECGP4uloU3qH~uCMH42MHf4TNnDL2zAsKA7cXsCsQYz0A9yJdNh7dfRT8ohpuAzASFx5Pj-bkqfw4p2tql55IIaPN0zqxyUuyZ9sfKl5qTQX82LoVolRpiBUL8dF9nr~bA2P1gJZ~xg1QssS7icR5MmTzvKKS5NYkezG8vWaTiEdXU0nuKI2ciZSX5GOMeIRW-YYR7mJwHmBbTVxe0o-uBuUtqor0Y98jdIv8Z~dwMjujRjrEShdCGNixTSonGzeS2~9CXqWquCJIOolqFFkcFHgXkD7ZWNfSXWbTxuF57rCsub98pA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

Reference: human and mouse
wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0.tar.gz

tar -xzvf cellranger-6.1.2.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0.tar.gz
```

## Usage

```
conda activate celescope
multi_flv_CR \
    --mapfile ./test.mapfile \
    --thread 8 \
    --seqtype TCR \
    --ref_path "/soft/cellranger/vdj_ref/6.0.0/hs/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0" \
    --soft_path "/soft/cellranger/cellranger-6.1.2/cellranger" \
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


### convert

- Convert barcodes and UMI to 10X format.

Output

- `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

- `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads as cellranger input.

- `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads as cellranger input.

### assemble

- TCR/BCR Assemble by Cellranger.

- Generate Mapping, Cells, V(D)J annotations metrics in html.


### summarize

- Convert 10X barcode of assemble result back to SGR barcode.

- Generate Productive contigs sequences and annotation files.

- Generate VDJ-annotation metrics in html.


### match

- Assembled results match with sc-RNA library.

- Generate matched VDJ-annotation metrics, clonetypes table and bar-plot of clonetypes distribution in html.


### mapping

- Output assembled T/B cells mapping to transcriptome if rds and auto-assign info exist in match directory.


## Output files
### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### assemble

- `03.assemble/{sample}` Cellranger vdj results.

### summarize

- `filtered_contig_annotations.csv` High-level annotations of each high-confidence contigs from cell-associated barcodes.

- `filtered_contig.fasta` High-confidence contig sequences annotated in the filtered_contig_annotations.csv.

- `productive_contig_annotations.csv` Annotations of each productive contigs from cell-associated barcodes. This is a subset of filtered_contig_annotations.csv.

- `productive_contig.fasta` Productive contig sequences annotated in the productive_contig_annotations.csv.

- `clonotypes.csv` High-level descriptions of each clonotype.

### match

- `matched_contig_annotations.csv` High-level annotations of each high-confidence contigs from matched cell-associated barcodes.

- `matched_contig.fasta` High-confidence contig sequences annotated in the matched_contig_annotations.csv.

- `matched_productive_contig_annotations.csv` Annotations of each productive contigs from matched cell-associated barcodes. This is a subset of matched_contig_annotations.csv.

- `matched_productive_contig.fasta` Productive contig sequences annotated in the matched_productive_contig_annotations.csv.

- `clonotypes.csv` High-level descriptions of each clonotype where barcodes match with scRNA-Seq.

### mapping
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

`--tenX_chemistry` 10X chemistry version, V2 or V3 for scRNA, V2 for VDJ.

`--ref_path` reference path for cellranger.

`--other_param` Other cellranger parameters.

`--mem` memory(G).

`--soft_path` soft path for cellranger.

`--seqtype` TCR or BCR.

