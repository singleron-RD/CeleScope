## Features
- Split scRNA-Seq fastq according to tag assignment.

## Output
- `matrix/` Matrix files of each tag.(Optional)
- `fastq/` Fastq files of each tag.(Optional)


## Arguments
`--split_fastq` If used, will split scRNA-Seq fastq file according to tag assignment.

`--split_matrix` If used, will split scRNA-Seq matrix file according to tag assignment.

`--umi_tag_file` UMI tag file.

`--match_dir` Match celescope scRNA-Seq directory.

`--matrix_dir` Match celescope scRNA-Seq matrix directory.

`--R1_read` R1 read path.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

