## Features
- Cell clustering with Seurat.

- Calculate the marker gene of each cluster.

- Cell type annotation(optional). You can provide markers of known cell types and annotate cell types for each cluster.

## Output
- `markers.tsv` Marker genes of each cluster.

- `tsne_coord.tsv` t-SNE coordinates and clustering information.

- `{sample}/06.analsis/{sample}_auto_assign/` This result will only be obtained when `--type_marker_tsv` 
parameter is provided. The result contains 3 files:
    - `{sample}_auto_cluster_type.tsv` The cell type of each cluster; if cell_type is "NA", 
it means that the given marker is not enough to identify the cluster.
    - `{sample}_png/{cluster}_pctdiff.png` Percentage of marker gene expression in this cluster - percentage in all other clusters.
    - `{sample}_png/{cluster}_logfc.png` log2 (average expression of marker gene in this cluster / average expression in all other clusters + 1)
## Arguments
`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--matrix_file` Required. Matrix_10X directory from step count.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

