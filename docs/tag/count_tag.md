## Features
- Assign tag to each cell barcode and summarize.

## Output

- `{sample}_umi_tag.tsv` 

    `first column` cell barcode  
    `last column`  assigned tag  
    `columns between first and last` UMI count for each tag 

- `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation

- `{sample}_cluster_count.tsv` cell barcode number assigned to *undeterminded*, *multiplet* and *each tag*
## Arguments
`--UMI_min` Default='auto'. Minimum UMI threshold. Cell barcodes with valid UMI < UMI_min are classified as *undeterminded*.

`--dim` Default=1. Tag dimentions. Usually we use 1-dimentional tag.

`--SNR_min` Default='auto'. Minimum signal-to-noise ratio. 
Cell barcodes with UMI >=UMI_min and SNR < SNR_min are classified as *multiplet*.

`--combine_cluster` Conbine cluster tsv file.

`--coefficient` Default=0.1. If `SNR_min` is 'auto', minimum signal-to-noise ratio is calulated as 
`SNR_min = max(median(SNRs) * coefficient, 2)`. 
Smaller `coefficient` will cause less *multiplet* in the tag assignment.

`--read_count_file` Tag read count file.

`--match_dir` Match celescope scRNA-Seq directory.

`--matrix_dir` Match celescope scRNA-Seq matrix directory.

`--tsne_file` match_dir t-SNE coord file. Do not required when `--match_dir` is provided.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

