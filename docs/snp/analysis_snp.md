## Features
- Annotate variants with [Annovar](https://annovar.openbioinformatics.org/en/latest/).

## Output
- `{sample}_gt.csv` Genotypes of variants of each cell. Rows are variants and columns are cells.
- `{sample}_variant_ncell.csv` Number of cells with each genotype.
- `{sample}_variant_table.csv` `{sample}_variant_ncell.csv` annotated with COSMIC(https://cancer.sanger.ac.uk/cosmic).
## Arguments
`--annovar_config` ANNOVAR config file.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` Match celescope scRNA-Seq directory.

`--vcf` vcf file.

