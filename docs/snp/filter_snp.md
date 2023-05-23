## Features
- Filter out `ref` and `alt` alleles that do not have enough reads to support.

## Output
- `{sample}_test1_filtered.vcf` VCF file after filtering. Alleles read counts that do not have enough reads to support are set to zero. 
Genotypes are changed accordingly.
## Arguments
`--ref_threshold_method` One of [otsu, auto, hard, none].

`--alt_threshold_method` One of [otsu, auto, hard, none].

`--ref_min_support_read` minimum supporting read number for ref.

`--alt_min_support_read` minimum supporting read number for alt.

`--vcf` norm vcf file.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

