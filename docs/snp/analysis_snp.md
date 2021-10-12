## Features
- Annotate variants with [Annovar](https://annovar.openbioinformatics.org/en/latest/).

## Output
- `{sample}.{genome_version}_multianno.txt` Annovar main output file. `CID` and `VID` are added to the `Otherinfo` column.

- `{sample}_variant_table.tsv` Formatted `multianno` file with `nCell`(number of cells with the variant) added.


## Arguments
`--annovar_config` ANNOVAR config file.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_dir` Match celescope scRNA-Seq directory.

`--filter_vcf` filter vcf file.

`--CID_file` CID_file.

`--filter_variant_count_file` filter variant count file.

`--ncell_file` filter cell count file.

