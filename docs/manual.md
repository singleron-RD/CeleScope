## Introduction
CeleScope is a collection of bioinfomatics analysis pipelines to process single cell sequencing data generated with Singleron products. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis as well as a summary of QC criteria.

Each pipeline consists of several steps and they all have two identical pre-processing steps: `barcode` and `cutadapt`. `barcode`step is used for barcode demupltiplexing, correction and read filtering. `cutadapt` step calls [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for read trimming.

Currently, CeleScope includes the follwing pipelines:

- `celescope rna` for Single-cell RNA-seq data generated with GEXSCOPE<sup>TM</sup> kits. It performs preprocessing, genome alignment, feature counting, expression matrix generation, clustering, marker gene expression analysis and cell type assignment(optional).

- `celescope vdj` for Single-cell Immune Repertoire data generated with GEXSCOPE<sup>TM</sup> IR kits. It performs preprocessing, UMI consensus, vdj sequence alignment, UMI filtering and clonetypes counting.

- `celescope tag` for Single-cell Multiplexing data generated with CLindex<sup>TM</sup> Sample Multiplexing kits. It performs preprocessing, tag counting, tag assignment and multiplets identification.

- `celescope dynaseq` for Single cell dynamic transcriptome data generated with DynaSCOPE<sup>TM</sup> kits. It performs preprocessing, genome alignment, feature counting, expression matrix generation, clustering, nucleotide substitution statistic and nascent RNA identification.

- `celescope snp` for Single cell variant data generated with FocuSCOPE<sup>TM</sup> kits. It performs preprocessing, genome alignment, feature counting, variant calling, variant filtering and variant annotation.

- `celescope capture_virus` for Single cell Epstein-Barr virus(EBV) data generated with FocuSCOPE<sup>TM</sup> Single Cell mRNA × EBV Library  kits. It performs preprocessing, genome alignment, UMI counting, filtering and analysis.

## [Installation](installation.md)

## [Quick start](quick_start.md)

## [Chemistry](chemistry.md)

## [Change log](CHANGELOG.md)

## [Test scripts and data](https://github.com/singleron-RD/celescope_test_script)

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
## Single-cell dynaseq
- [star](rna/star.md)
- [featureCounts](tools/featureCounts.md)
- [count](tools/count.md)
- [analysis](rna/analysis.md)
- [conversion](dynaseq/conversion.md)
- [substitution](dynaseq/substitution.md)
- [replacement](dynaseq/replacement.md)
- [replace_tsne](dynaseq/replace_tsne.md)
- [multi_dynaseq](dynaseq/multi_dynaseq.md)
## Single-cell snp
- [mkref](snp/mkref.md)
- [consensus](tools/consensus.md)
- [star](rna/star.md)
- [featureCounts](tools/featureCounts.md)
- [target_metrics](tools/target_metrics.md)
- [variant_calling](snp/variant_calling.md)
- [filter_snp](snp/filter_snp.md)
- [analysis_snp](snp/analysis_snp.md)
- [multi_snp](snp/multi_snp.md)
## Single-cell capture_virus
- [mkref](capture_virus/mkref.md)
- [consensus](tools/consensus.md)
- [star_virus](rna_virus/star_virus.md)
- [count_virus](capture_virus/count_virus.md)
- [filter_virus](capture_virus/filter_virus.md)
- [analysis_virus](capture_virus/analysis_virus.md)
- [multi_capture_virus](capture_virus/multi_capture_virus.md)
