## Features
- Annotate variants with [snpEff](http://pcingola.github.io/SnpEff/).

## Output
- `{sample}_gt.csv` Genotypes of variants of each cell. Rows are variants and columns are cells.
- `{sample}_variant_ncell.csv` Number of cells with each genotype.
- `{sample}_variant_table.csv` annotated with snpEff.
## Arguments
`--gene_list` Required. Gene list file, one gene symbol per line. Only results of these genes are reported. Conflict with `--panel`.

`--database` snpEff database. Common choices are GRCh38.99(human) and GRCm38.99(mouse).

`--panel` The prefix of bed file in `celescope/data/snp/panel/`, such as `lung_1`. Conflict with `--gene_list`.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` Match celescope scRNA-Seq directory.

`--vcf` vcf file.

