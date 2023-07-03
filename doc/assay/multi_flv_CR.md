## Download and unpack cellranger soft and reference file.
```
wget -O cellranger-7.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-vdj/cellranger-7.0.1.tar.gz?Expires=1668447380&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC12ZGovY2VsbHJhbmdlci03LjAuMS50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2Njg0NDczODB9fX1dfQ__&Signature=bZoM-0W1ExsZWYVhUJ80oSzKEaYTdlf1NoFiKLswTLWGveefjH1fJOkd1c4kMxCiTDfohyDpNWzk-xBMBme-u3r-0X7WFYSvfCdMyo2CSM4K~Ur73Sn30REr0kr7oC9byrlLMqG2mtO7DLfprDrAZGZDkfyDLpGJ8hb1qWSsWqRo8CxbqRzA69h0v65Qn86HMrHAeotdhxUVmBIvONwPmsC90J9K4gVDw1sDF39F4f89zguTGgSAY6mUdPG1cvHHZMeuLaJZDqRgODysOhsB-keYXW8cYa-R5chh9s1ASC4yA1QKSe-fBdB1-FLQdAwMjNEHdd7uGCWzLj4J~AXgJg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

Reference: human and mouse
wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz

tar -xzvf cellranger-7.0.1.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0.tar.gz
tar -xzvf refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz
```

## Usage

```
conda activate celescope
multi_flv_CR \
    --mapfile ./test.mapfile \
    --thread 8 \
    --seqtype TCR \
    --ref_path "/soft/cellranger/vdj_ref/6.0.0/hs/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0" \
    --soft_path "/soft/cellranger/cellranger-7.0.1/cellranger" \
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

- Output TSNE-plot of Assembled T/B Cells.


### refine_vdj

- Refine barcodes where "is_cell=False" and have multi productive chains.

- There are three methods to filter noise: SNR, AUTO, NOT_FILTER

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
- `06.mapping/{sample}_mapping.pdf` TSNE-plot of Assembled Cells.

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

`--tenX_chemistry` 10X chemistry version, V2 or V3 for scRNA, V2 for VDJ.

`--ref_path` reference path for cellranger.

`--other_param` Other cellranger parameters.

`--mem` memory(G).

`--soft_path` soft path for cellranger.

`--not_refine` Run the VDJ refine step.

`--filter_method` filter noise method.

`--coeff` coefficient will affect auto and snr noise filter, recommend 1.5 for auto, 10 for snr.

`--seqtype` TCR or BCR.

