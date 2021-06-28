## Features
- Split scRNA-Seq fastq according to tag assignment.

## Output
- `fastq/{tag}_{1,2}.fq` Fastq files of each tag.


## Arguments
`--split_fastq` If used, will split scRNA-Seq fastq file according to tag assignment.

`--umi_tag_file` UMI tag file.

`--match_dir` Match celescope scRNA-Seq directory.

`--R1_read` R1 read path.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

