## Features
- Perform variant calling at single cell level.

## Output

- `{sample}_norm.vcf` Normalized vcf file.


## Arguments
`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--panel` The prefix of bed file in `celescope/data/snp/panel/`, such as `lung_1`.

`--bam` Input BAM file from step `target_metrics`.

`--match_dir` Match celescope scRNA-Seq directory.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

