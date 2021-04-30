# STAR

## Features
- Align R2 reads to the reference genome with STAR.
- Collect Metrics with Picard.

## Input
- R2 clean reads from step cutadapt.

## Output
- `{sample}_Aligned.sortedByCoord.out.bam` BAM file contains Uniquely Mapped Reads.

- `{sample}_SJ.out.tab` SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format.

- `{sample}_Log.out` Main log with a lot of detailed information about the run. This is most useful for troubleshooting and debugging.

- `{sample}_Log.progress.out` Report job progress statistics, such as the number of processed reads, % of mapped reads etc. It is updated in 1 minute intervals.

- `{sample}_Log.Log.final.out` Summary mapping statistics after mapping job is complete, very useful for quality control. The statistics are calculated for each read (single- or paired-end) and then summed or averaged over all reads. Note that STAR counts a paired-end read as one read, (unlike the samtools agstat/idxstats, which count each mate separately). Most of the information is collected about the UNIQUE mappers (unlike samtools agstat/idxstats which does not separate unique or multi-mappers). Each splicing is counted in the numbers of splices, which would correspond to summing the counts in SJ.out.tab. The mismatch/indel error rates are calculated on a per base basis, i.e. as total number of mismatches/indels in all unique mappers divided by the total number of mapped bases.

- `{sample}_region.log` Picard CollectRnaSeqMetrics results.

## Parameters

`--fq` Required. R2 reads from step cutadapt.

`--genomeDir` Directory contains genome Fasta, GTF and refFLAT file. If this argument is not provided, you need to provide `--STAR_index` and `--refFlat`

`--STAR_index` STAR index directory path.

`--refFlat` refFlat file path.

`--starMem` Default `30`. Maximum memory that STAR can use.

`--out_unmapped` will output unmapped and partially mapped (i.e. mapped only one
mate of a paired end read) reads into separate Unmapped.out.mate1(2), formatted the same
way as input read (i.e. FASTQ or FASTA).

`--outFilterMultimapNmax` Default `1`. How many places are allowed to match a read at most.

`--outFilterMatchNmin` Default `0`. Alignment will be output only if the number of matched bases is higher than or equal to this value.

`--STAR_param` Default `''`. Other STAR parameters.

## Metrics
- Uniquely Mapped Reads: reads that mapped uniquely to the genome.

- Multi-Mapped Reads: reads that mapped to multiple locations in the genome.

- Base Pairs Mapped to Exonic Regions: bases in primary alignments that align to a coding base or a UTR base for some gene.

- Base Pairs Mapped to Intronic Regions: bases in primary alignments that align to an intronic base for some gene, and not a coding or UTR base.

- Base Pairs Mapped to Intergenic Regions: bases in primary alignments that do not align to any gene.