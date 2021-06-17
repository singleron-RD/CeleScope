## Features
- Perform variant calling.

## Output

`{sample}_VID.tsv` A unique numeric ID is assigned for each variant.

`{sample}_CID.tsv` A unique numeric ID is assigned for each cell.

`{sample}_variant_count.tsv`  Reference and variant supporting reads/UMIs count.

`{sample}_support.mtx` Support matrix, only high quality bases are considered.   
0 : no reads/UMIs cover the position.  
1 : all reads/UMIs at the position support the ref allele.  
2 : all reads/UMIs at the position support the alt allele.  
3 : one or more reads/UMIs support both the alt and the ref allele.  


## Arguments
`--genomeDir` Genome directory after running `mkref`.

`--vcf` VCF file. If vcf file is not provided, celescope will perform variant calling at single cell level 
and use these variants as input vcf.

`--bam` Input BAM file from step `target_metrics`.

`--match_dir` Match celescope scRNA-Seq directory.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

