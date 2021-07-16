## Features
- Filter bam file
    - Filter reads that are not cell-associated.
    - Filter reads that are not mapped to target genes. 

- Collect enrichment metrics.

## Output
- `filtered.bam` BAM file after filtering.


## Arguments
`--gene_list` Required. Gene list file, one gene symbol per line. Only results of these genes are reported.

`--bam` Input bam file

`--match_dir` Match celescope scRNA-Seq directory.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

