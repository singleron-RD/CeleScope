## Features
- Perform variant calling at single cell level.

## Output
- `{sample}_raw.vcf` Variants are called with bcftools default settings.
- `{sample}_norm.vcf` Indels are left-aligned and normalized. See https://samtools.github.io/bcftools/bcftools.html#norm for more details.
## Arguments
`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--panel` The prefix of bed file in `celescope/data/snp/panel/`, such as `lung_1`. Conflict with `--gene_list`.

`--bam` Input BAM file from step `target_metrics`.

`--match_dir` Match celescope scRNA-Seq directory.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

