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
`--genomeDir` Required. Genome directory.

`--save_rds` Write rds to disk.

`--type_marker_tsv` A tsv file with header. If this parameter is provided, cell type will be annotated. Example:
```
cell_type	marker
Alveolar	"CLDN18,FOLR1,AQP4,PEBP4"
Endothelial	"CLDN5,FLT1,CDH5,RAMP2"
Epithelial	"CAPS,TMEM190,PIFO,SNTN"
Fibroblast	"COL1A1,DCN,COL1A2,C1R"
B_cell	"CD79A,IGKC,IGLC3,IGHG3"
Myeloid	"LYZ,MARCO,FCGR3A"
T_cell	"CD3D,TRBC1,TRBC2,TRAC"
LUAD	"NKX2-1,NAPSA,EPCAM"
LUSC	"TP63,KRT5,KRT6A,KRT6B,EPCAM"
```

`--matrix_file` Required. Matrix_10X directory from step count.

`--tsne_file` match_dir t-SNE coord file.

`--df_marker_file` match_dir df_marker_file.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

