## Features
- Assign tag to each cell barcode and summarize.

## Output

- `{sample}_umi_tag.tsv` 

    `first column` cell barcode  
    `last column`  assigned tag  
    `columns between first and last` UMI count for each tag 

- `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation
## Arguments
`--read_count_file` Tag read count file.

`--match_dir` Match celescope scRNA-Seq directory.

`--matrix_dir` Match celescope scRNA-Seq matrix directory.

`--tsne_file` match_dir t-SNE coord file. Do not required when `--match_dir` is provided.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

