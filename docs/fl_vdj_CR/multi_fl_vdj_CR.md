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
multi_fl_vdj_CR \
    --mapfile ./test.mapfile \
    --chemistry flv \
    --mem 10 \
    --thread 8 \
    --allowNoLinker \
    --seqtype TCR \
    --ref_path "/SGRNJ/Database/script/soft/cellranger/vdj_ref/6.0.0/hs/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0" \
    --soft_path "/SGRNJ03/randd/cjj/soft/cellranger/cellranger-6.1.2/cellranger" 
```
## Features
### barcode

- Demultiplex barcodes.
- Filter invalid R1 reads, which includes:
    - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
    - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
    - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
    - Low quality reads: low sequencing quality in barcode and UMI regions.


### assemble

- TCR/BCR Assemble.


### annotation

- V(D)J annotation infomation.


### match

- V(D)J results match SC-RNA infomation.


### mapping

- Assembled T/B cells Mapping with SC-RNA barcodes.


## Output files
### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### assemble
- `03.assemble/{sample}/outs/` Recording assemble results.

### annotation
- `filtered_contig_annotations.csv' High-level annotations of each high-confidence, cellular contig.

- `filtered_contig.fasta` High-confidence contig sequences in cell barcodes.

- `clonotypes.csv` High-level descriptions of each clonotype.

### match
- `match_contigs.csv` Consider barcodes match scRNA-Seq library in filtered_contig_annotations.csv.

- `match_contig.fasta` Consider barcodes match scRNA-Seq library in filtered_contig.fasta.

- `match_clonotypes.csv` Consider barcodes match scRNA-Seq library in clonotypes.csv.

### mapping
- `{sample}_assign.png` Auto-assigned umap plot in scRNA-Seq library.

- `{sample}_cluster_umap.png` Cluster umap plot in scRNA-Seq library.

- `{sample}_umapplot.png` Umap plot of assembled and un-assembled barcodes in scRNA-Seq library.

- `{sample}_distribution.txt` Number of assembled barcodes in every clusters in scRNA-Seq library.

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

`--gzip` Output gzipped fastq files.

`--output_R1` Output valid R1 reads.

`--species` species.

`--soft` cellranger version.

`--ref_path` reference path for cellranger.

`--soft_path` soft path for cellranger.

`--other_param` Other cellranger parameters.

`--not_split_R2` whether split r2.

`--mem` memory (G).

`--seqtype` TCR or BCR.

