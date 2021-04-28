# featureCounts

## Features
- Quantify genes in BAM files.

## Input
- BAM file from step STAR.

## Output
- `{sample}` Numbers of reads assigned to features (or meta-features).

- `{sample}_summary` Stat info for the overall summrization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons (these reasons are included in the stat info).

- `{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam` featureCounts output BAM, sorted by coordinates；BAM file contains tags as following(Software Version>=1.1.8):
CB cell barcode
UB UMI
GN gene name
GX gene id

- `{sample}_name_sorted.bam` featureCounts output BAM, sorted by read name.

## Parameters

`--genomeDir` Required. Directory contains genome Fasta, GTF and refFLAT file.

`--gtf_type` Default "exon". Single-nuclei RNA-seq uses 'gene'.

`--input` BAM file path.

## Metrics
- Assigned : reads that can be successfully assigned without ambiguity.

- Unassigned_NoFeatures : alignments that do not overlap any feature.

- Unassigned_Ambiguity : alignments that overlap two or more features.

