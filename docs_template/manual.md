## Introduction
CeleScope is a collection of bioinfomatics analysis pipelines developed at Singleron to process single cell sequencing data generated with Singleron products. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis as well as a summary of QC criteria.

Each pipeline consists of several steps and they all have two identical pre-processing steps: `barcode` and `cutadapt`. `barcode`step is used for barcode demupltiplexing, correction and read filtering. `cutadapt`step calls [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for read trimming.

Currently, CeleScope includes the follwing pipelines:

- `celescope rna` for Single-cell RNA-seq data generated with GEXSCOPE kits. It performs preprocessing, genome alignment, feature counting, expression matrix generation, clustering, marker gene expression analysis and cell type assignment(optional).

- `celescope vdj` for Single-cell Immune Repertoire data generated with GEXSCOPE IR kits. It performs preprocessing, UMI consensus, vdj sequence alignment, UMI filtering and clonetypes counting.

- `celescope tag` for Single-cell Multiplexing data generated with CLindex Sample Multiplexing kits. It performs preprocessing, tag counting, tag assignment and multiplets identification.

- `celescope dynaseq` for Single cell dynamic transcriptome data generated with Dynascope kits. It performs preprocessing, genome alignment, feature counting, expression matrix generation, clustering, nucleotide substitution statistic and nascent RNA identification.


## [Quick start](quick_start.md)

## [Change log](CHANGELOG.md)

## Pre-processing

- [barcode](tools/barcode.md)
- [cutadapt](tools/cutadapt.md)

