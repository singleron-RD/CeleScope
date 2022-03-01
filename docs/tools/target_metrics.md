## Features
- Filter bam file
    - Filter reads that are not cell-associated.
    - Filter reads that are not mapped to target genes. 

- Collect enrichment metrics.

## Output
- `filtered.bam` BAM file after filtering. Keep reads that are cell-associated and mapped to target genes.
## Arguments
`--gene_list` Required. Gene list file, one gene symbol per line. Only results of these genes are reported.

`--panel` The prefix of bed file in `celescope/data/snp/panel/`, such as `lung_1`.

`--bam` Input bam file.

`--match_dir` Match celescope scRNA-Seq directory.

`--add_RG` Add tag read group: RG. RG is the same as CB(cell barcode).

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

