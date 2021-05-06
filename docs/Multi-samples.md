# Multi-samples and step-by-step

## Features
- CeleScope contains an interface to run multiple samples, which is `multi_{assay}`.
- Step-by-step analysis.

## Input
- Mapfile. Mapfile is a text file divided by tabs, containing three columns, each line represents a sample.
1st column: Fastq prefix.
2nd column: Folder where fastq is located.
3rd column: {sample} (sample name, which is the prefix of all generated files).
4th column: optional, expected cell number (scRNA-Seq) or match_dir (scVDJ).
Note: When a sample has multiple fastqs, each fastq occupies one line, and the sample name is the same. During processing, these fastqs are treated as the same sample and merged.
	- Example:
	```
	$cat ./my.mapfile
	R2007197	/SGRNJ/DATA_PROJ/2003/20200728	S20071508_D_TS
	R2007198	/SGRNJ/DATA_PROJ/2003/20200728	S20071509_D_TS
	R2007197        /SGRNJ/DATA_PROJ/2003/test      S20071508_D_TS

	$ls -l /SGRNJ/DATA_PROJ/2003/20200728
	-rw-r--r--. 1 download ssh.bioinfo 1.2G Jul 29 00:50 R2007198_L2_2.fq.gz
	-rw-r--r--. 1 download ssh.bioinfo 1.3G Jul 29 00:21 R2007198_L2_1.fq.gz
	-rw-r--r--. 1 download ssh.bioinfo 1.3G Jul 28 23:52 R2007197_L2_2.fq.gz
	-rw-r--r--. 1 download ssh.bioinfo 1.4G Jul 28 23:23 R2007197_L2_1.fq.gz
	```

## Parameters

`--mod` Default `sjm`. Out file mode, `sjm` or `shell`. 

`--mapfile` Required. Mapfile path.

`--rm_files` Remove redundant fq.gz and bam after running.

`--steps_run` Default `all`. Steps to run.

`--outdir` Default `./`. Output dir.

`--debug` Debug or not.

`--thread`  Default `4`. Thread to use.

## Usage
- Multi-sample example(Single Cell RNA-Seq):

	```
	conda activate celescope
	multi_rna\
 	--mapfile ./my.mapfile\
 	--chemistry auto\
 	--genomeDir /SGRNJ/Public/Database/genome/homo_mus\
 	--thread 8\
 	--mod shell
 	```

 - Note: The recommended setting for thread is 8, and the maximum should not exceed 20.

 Scripts above will generate `shell` directory containing `{sample}.sh` file.

 You can start your analysis by running:
 ```
 sh ./shell/S20071508_D_TS.sh
 ```

 - Step-by-step:
 `{sample}.sh` records analysis scripts for each step. You can run these steps respectively. For example:
 ```
 celescope rna STAR --fq .//S20071508_D_TS/02.cutadapt/S20071508_D_TS_clean_2.fq.gz --sample S20071508_D_TS --genomeDir /SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92 --thread 6 --outdir .//S20071508_D_TS/03.STAR --assay rna
 ```
 
