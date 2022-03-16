## Features
- Computes the replacement rates in each cell and gene.
- Boxplots for rates distribution.

## Output
- `{sample}.TC_matrix.rds` New and old info for each barcode/gene/umi.
- `{sample}.new_matrix.tsv.gz` New RNA matrix.
- `{sample}.old_matrix.tsv.gz` Old RNA matrix.
- `{sample}.fraction_of_newRNA_per_cell.txt` Fraction of new RNA of each cell.
- `{sample}.fraction_of_newRNA_per_gene.txt` Fraction of new RNA of each gene.
- `{sample}.fraction_of_newRNA_matrix.txt` Fraction of new RNA of each cell and gene.
## Arguments
`--bg_cov` background snp depth filter, lower than bg_cov will be discarded. Only valid in csv format.

`--bam` bam file from conversion step.

`--bg` background snp file, csv or vcf format.

`--cell_keep` filter cell.

`--min_cell` a gene expressed in at least cells, default 10.

`--min_gene` at least gene num in a cell, default 10.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

