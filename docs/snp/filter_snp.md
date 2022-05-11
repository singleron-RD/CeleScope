## Features
- Filter out `ref` and `alt` alleles that do not have enough reads to support.

## Output
- `{sample}_test1_filtered.vcf` VCF file after filtering. Alleles read counts that do not have enough reads to support are set to zero. 
Genotypes are changed accordingly.
## Arguments
`--threshold_method` One of [otsu, auto, hard, none].

`--hard_threshold` int, use together with `--threshold_method hard`.

`--vcf` norm vcf file.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

