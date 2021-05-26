CeleScope is a collection of bioinfomatics analysis pipelines developed at Singleron to process single cell sequencing data generated with Singleron products. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis as well as a summary of QC criteria.

Each pipeline consists of several steps and they all have two identical preprocessing steps: `barcode` and `cutadapt`. `barcode`step is used for barcode demupltiplexing, correction and read filtering. `cutadapt`step calls [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for read trimming.

Currently, CeleScope includes the follwing pipelines:

- `celescope rna` for Single-cell RNA-seq data generated with GEXSCOPE kits. It performs preprocessing, genome alignment, feature counting, expression matrix generation, clustering, marker gene expression analysis and cell type assignment(optional).

- `celescope vdj` for Single-cell Immune Repertoire data generated with GEXSCOPE IR kits. It performs preprocessing, UMI consensus, vdj sequence alignment, UMI filtering and clonetypes counting. 

- `celescope tag` for Single-cell Multiplexing data generated with CLindex Sample Multiplexing kits. It performs preprocessing, tag counting, tag assignment and multiplets identification.

Celescope has interfaces `multi_{assay}` to generate scripts for multi-sample analysis. See `Multi-samples` for details.

- Common Steps
    - [barcode](tools/barcode.md)
    - [cutadapt](tools/cutadapt.md)

- Single cell RNA-Seq
    - [mkref](rna/mkref.md)
    - [star](rna/star.md)
    - [featureCounts](tools/featureCounts.md)
    - [count](rna/count.md)
    - [analysis](rna/analysis.md)

- Single cell VDJ
    - [consensus](tools/consensus.md)
    - [mapping_vdj](vdj/mapping_vdj.md)
    - [count_vdj](vdj/count_vdj.md)

- [Multi-samples](multi-samples.md)

- [Change log](CHANGELOG.md)