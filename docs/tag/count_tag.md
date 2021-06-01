# count_tag

## Features
- Assign tag to each cell barcode and summarize

## Input
- Read count file from mapping_tag
- matched rna directory

## Output

- `{sample}_umi_tag.tsv` 

    `first column` cell barcode  
    `last column`  assigned tag  
    `columns between first and last` UMI count for each tag 

- `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation

- `{sample}_cluster_count.tsv` cell barcode number assigned to *undeterminded*, *multiplet* and *each tag*


## Parameters

`--UMI_min` Default='auto'. Minimum UMI threshold. Cell barcodes with valid UMI < UMI_min are classified as *undeterminded*.

`--SNR_min` Default='auto'. Minimum signal-to-noise ratio. Cell barcodes with UMI >=UMI_min and SNR < SNR_min are classified as *multiplet*.

`--dim` Default=1. Tag dimentions. Usually we use 1-dimentional tag.

`--coefficient` Default=0.1. If `SNR_min` is 'auto', minimum signal-to-noise ratio is calulated as `SNR_min = max(median(SNRs) * coefficient, 2)`. Smaller `coefficient` will cause less *multiplet* in the tag assignment.

## Metrics

- Mapped Reads in Cells : Mapped reads with scRNA-Seq cell barcode

- Median UMI per Cell : Median UMI per scRNA-Seq cell barcode

- Mean UMI per Cell : Mean UMI per scRNA-Seq cell barcode