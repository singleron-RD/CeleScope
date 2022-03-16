## Features
- Align R2 reads to the reference genome with STAR.
- Collect Metrics with Picard.

## Output
- `{sample}_Aligned.sortedByCoord.out.bam` BAM file contains Uniquely Mapped Reads.

- `{sample}_SJ.out.tab` SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format.

- `{sample}_Log.out` Main log with a lot of detailed information about the run. 
This is most useful for troubleshooting and debugging.

- `{sample}_Log.progress.out` Report job progress statistics, such as the number of processed reads, 
% of mapped reads etc. It is updated in 1 minute intervals.

- `{sample}_Log.Log.final.out` Summary mapping statistics after mapping job is complete, 
very useful for quality control. The statistics are calculated for each read (single- or paired-end) and 
then summed or averaged over all reads. Note that STAR counts a paired-end read as one read, 
(unlike the samtools agstat/idxstats, which count each mate separately). 
Most of the information is collected about the UNIQUE mappers 
(unlike samtools agstat/idxstats which does not separate unique or multi-mappers). 
Each splicing is counted in the numbers of splices, which would correspond to 
summing the counts in SJ.out.tab. The mismatch/indel error rates are calculated on a per base basis, 
i.e. as total number of mismatches/indels in all unique mappers divided by the total number of mapped bases.

- `{sample}_region.log` Picard CollectRnaSeqMetrics results.
## Arguments
`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--outFilterMatchNmin` Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.

`--out_unmapped` Output unmapped reads.

`--STAR_param` Other STAR parameters.

`--outFilterMultimapNmax` Default `1`. How many places are allowed to match a read at most.

`--starMem` Default `30`. Maximum memory that STAR can use.

`--fq` Required. R2 fastq file.

`--consensus_fq` A indicator that the input fastq has been consensused.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

