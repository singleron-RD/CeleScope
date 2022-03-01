## Features
- Annotate variants with [Annovar](https://annovar.openbioinformatics.org/en/latest/).

## Output
- `{sample}_gt.csv` genotypes of variants of each cell. Row is variant, column is cell.
- `{sample}_variant_ncell.csv` Number of cells with each genotype.
- `{sample}_variant_table.csv` Annotated `{sample}_variant_ncell.csv`.
## Arguments
`--annovar_config` ANNOVAR config file.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` Match celescope scRNA-Seq directory.

`--vcf` vcf file.

