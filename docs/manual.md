## Introduction
CeleScope is a collection of bioinfomatics analysis pipelines developed at Singleron to process single cell sequencing data generated with Singleron products. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis as well as a summary of QC criteria.

Each pipeline consists of several steps and they all have two identical pre-processing steps: `barcode` and `cutadapt`. `barcode`step is used for barcode demupltiplexing, correction and read filtering. `cutadapt`step calls [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for read trimming.

Currently, CeleScope includes the follwing pipelines:

- `celescope rna` for Single-cell RNA-seq data generated with GEXSCOPE kits. It performs preprocessing, genome alignment, feature counting, expression matrix generation, clustering, marker gene expression analysis and cell type assignment(optional).

- `celescope vdj` for Single-cell Immune Repertoire data generated with GEXSCOPE IR kits. It performs preprocessing, UMI consensus, vdj sequence alignment, UMI filtering and clonetypes counting.

- `celescope tag` for Single-cell Multiplexing data generated with CLindex Sample Multiplexing kits. It performs preprocessing, tag counting, tag assignment and multiplets identification.


## [Quick start](quick_start.md)

## [Change log](CHANGELOG.md)

## Pre-processing

- [barcode](tools/barcode.md)
- [cutadapt](tools/cutadapt.md)

## Single-cell rna
- [mkref](rna/mkref.md)
- [star](rna/star.md)
- [featureCounts](tools/featureCounts.md)
- [count](tools/count.md)
- [analysis](rna/analysis.md)
- [multi_rna](rna/multi_rna.md)
## Single-cell vdj
- [consensus](tools/consensus.md)
- [mapping_vdj](vdj/mapping_vdj.md)
- [count_vdj](vdj/count_vdj.md)
- [multi_vdj](vdj/multi_vdj.md)
## Single-cell tag
- [mapping_tag](tag/mapping_tag.md)
- [count_tag](tag/count_tag.md)
- [analysis_tag](tag/analysis_tag.md)
- [split_tag](tag/split_tag.md)
- [multi_tag](tag/multi_tag.md)
## Single Cell Dynaseq
- [star](rna/star.md)
- [featureCounts](tools/featureCounts.md)
- [count](tools/count.md)
- [analysis](rna/analysis.md)
- [conversion](dynaseq/conversion.md)
- [subsitution](dynaseq/subsitution.md)
- [replacement](dynaseq/replacement.md)
- [replace_tsne](dynaseq/replace_tsne.md)
- [multi_dynaseq](dynaseq/multi_dynaseq.md)
