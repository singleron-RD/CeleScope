## Features
- Split scRNA-Seq fastq according to tag assignment.

## Output
- `matrix/` Matrix files of each tag.(Optional)
- `fastq/` Fastq files of each tag.(Optional)
## Arguments
`--split_fastq` If used, will split scRNA-Seq fastq file according to tag assignment.

`--split_matrix` If used, will split scRNA-Seq matrix file according to tag assignment.

`--split_bam` If used, will split scRNA-Seq featureCounts bam file according to tag assignment. Use together with --bam_file.

`--bam_file` scRNA-Seq bam file to split.

`--split_vdj` If used, will split scRNA-Seq vdj count file according to tag assignment.

`--split_fl_vdj` If used, will split scRNA-Seq full-length vdj annotation, fasta, clonotypes file according to tag assignment.

`--vdj_dir` Match celescope vdj directory. Required when --split_vdj or --split_fl_vdj is specified.

`--umi_tag_file` UMI tag file.

`--match_dir` Match celescope scRNA-Seq directory.

`--matrix_dir` Match celescope scRNA-Seq matrix directory.

`--R1_read` R1 read path.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

