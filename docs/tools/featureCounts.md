## Features
- Assigning uniquely mapped reads to genomic features with FeatureCounts.
## Output
- `{sample}` Numbers of reads assigned to features (or meta-features).
- `{sample}_summary` Stat info for the overall summrization results, including number of 
successfully assigned reads and number of reads that failed to be assigned due to 
various reasons (these reasons are included in the stat info).
- `{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam` featureCounts output BAM, 
sorted by coordinates；BAM file contains tags as following(Software Version>=1.1.8):
    - CB cell barcode
    - UB UMI
    - GN gene name
    - GX gene id
- `{sample}_name_sorted.bam` featureCounts output BAM, sorted by read name.
## Arguments
`--gtf_type` Specify feature type in GTF annotation.

`--genomeDir` Required. Genome directory.

`--featureCounts_param` Other featureCounts parameters.

`--input` Required. BAM file path.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

