## Features

## Output

- `{sample}_umi_tag.tsv` 

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

