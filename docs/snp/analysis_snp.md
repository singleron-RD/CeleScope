## Features
- Annotate variants with [Annovar](https://annovar.openbioinformatics.org/en/latest/).

## Output
`{sample}.{genome_version}_multianno.txt` Annovar main output file. `CID` and `VID` are added to the `Otherinfo` column.

`{sample}_variant_table.tsv` Formatted `multianno` file with `ncell_cover`and `ncell_alt` added.

`{sample}_variant_top5.jpg` The Venn diagram of the 5 variants with the highest `ncell_alt`.

`{sample}_variant_ncell.tsv` Number of cells with read count at each variant's position. 
- `VID`: Variant ID. 
- `ncell_cover`: number of cells with read count at this position. 
- `ncell_alt`: number of cells with variant read count only. 
- `ncell_ref`: number of cells with reference read count only. 
- `ncell_ref_and_alt`: number of cells with both variant and reference read count.
- `RID`: Target region ID. This column will be added when `--panel` option were provided.


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

