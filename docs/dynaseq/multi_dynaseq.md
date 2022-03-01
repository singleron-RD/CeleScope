## Usage
```
    multi_dynaseq\
    --mapfile ./rna.mapfile\
    --genomeDir /SGRNJ/Public/Database/genome/homo_mus\
    --STAR_param "--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outSAMattributes MD"\
    --strand /SGRNJ03/Public/Database/genome/gene.strandedness.csv
```

You need to generate strandness-file from gtf file. 
The format is "geneID,strand", eg:
```
ENSG00000223972,+
ENSG00000227232,-
ENSG00000278267,-
```
## Features
### barcode

- Demultiplex barcodes.
- Filter invalid R1 reads, which includes:
    - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
    - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
    - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
    - Low quality reads: low sequencing quality in barcode and UMI regions.


### cutadapt
- Trim adapters in R2 reads with cutadapt. Default adapters includes:
    - polyT=A{18}, 18 A bases. 
    - p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Illumina p5 adapter.

### star
- Align R2 reads to the reference genome with STAR.
- Collect Metrics with Picard.


### featureCounts
- Assigning uniquely mapped reads to genomic features with FeatureCounts.

### count
- Cell-calling: Distinguish cell barcodes from background barcodes. 
- Generate expression matrix.

### analysis
- Cell clustering with Seurat.

- Calculate the marker gene of each cluster.

- Cell type annotation(optional). You can provide markers of known cell types and annotate cell types for each cluster.


### conversion
- Get conversion pos in each read.
    - Get snp info. 


### substitution
- Computes the overall conversion rates in reads and plots a barplot.


### replacement
- Computes the replacement rates in each cell and gene.
- Boxplots for rates distribution.


### replace_tsne
- Replace rate in each cluster
- Top replace genes in each cluster


## Output files
### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
the read name is `{barcode}_{UMI}_{read ID}`.

### cutadapt
- `cutadapt.log` Cutadapt output log file.
- `{sample}_clean_2.fq.gz` R2 reads file without adapters.

### star
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

### featureCounts
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

### count
- `{sample}_all_matrix` The expression matrix of all detected barcodes. 
    Can be read in by calling the `Seurat::Read10X` function.
- `{sample}_matrix_10X` The expression matrix of the barcode that is identified to be the cell. 
Can be read in by calling the `Seurat::Read10X` function.
- `{sample}_matrix.tsv.gz` The expression matrix of the barcode that is identified to be the cell, separated by tabs. 
CeleScope >=1.2.0 does not output this file.
- `{sample}_count_detail.txt.gz` 4 columns: 
    - barcode  
    - gene ID  
    - UMI count  
    - read_count  
- `{sample}_counts.txt` 6 columns:
    - Barcode: barcode sequence
    - readcount: read count of each barcode
    - UMI2: UMI count (with reads per UMI >= 2) for each barcode
    - UMI: UMI count for each barcode
    - geneID: gene count for each barcode
    - mark: cell barcode or backgound barcode.
        `CB` cell  
        `UB` background  
- `{sample}_downsample.txt` 3 columns：
    - percent: percentage of sampled reads
    - median_geneNum: median gene number per cell
    - saturation: sequencing saturation
- `barcode_filter_magnitude.pdf` Barcode-UMI plot.

### analysis
- `markers.tsv` Marker genes of each cluster.

- `tsne_coord.tsv` t-SNE coordinates and clustering information.

- `{sample}/06.analsis/{sample}_auto_assign/` This result will only be obtained when `--type_marker_tsv` 
parameter is provided. The result contains 3 files:
    - `{sample}_auto_cluster_type.tsv` The cell type of each cluster; if cell_type is "NA", 
it means that the given marker is not enough to identify the cluster.
    - `{sample}_png/{cluster}_pctdiff.png` Percentage of marker gene expression in this cluster - percentage in all other clusters.
    - `{sample}_png/{cluster}_logfc.png` log2 (average expression of marker gene in this cluster / average expression in all other clusters + 1)

### conversion
- `{sample}.PosTag.bam` Bam file with conversion info.
- `{sample}.PosTag.csv` SNP info in csv format.

### substitution
- `{sample}.substitution.txt` Tab-separated table of the overall conversion rates.

### replacement
- `{sample}.TC_matrix.rds` New and old info for each barcode/gene/umi.
- `{sample}.new_matrix.tsv.gz` New RNA matrix.
- `{sample}.old_matrix.tsv.gz` Old RNA matrix.
- `{sample}.fraction_of_newRNA_per_cell.txt` Fraction of new RNA of each cell.
- `{sample}.fraction_of_newRNA_per_gene.txt` Fraction of new RNA of each gene.
- `{sample}.fraction_of_newRNA_matrix.txt` Fraction of new RNA of each cell and gene.

### replace_tsne
- `{sample}.rep_in_tsne.txt` Replace rate in each cluster.
- `{sample}.rep_in_tsne_top10` Top 10 replace genes in each cluster.

## Arguments
`--mapfile` Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  
4th column: The 4th column has different meaning for each assay. The single cell rna directory after running CeleScope is called `matched_dir`.

- `rna` Optional, forced cell number.
- `vdj` Optional, matched_dir.
- `tag` Required, matched_dir.
- `dynaseq` Optional, forced cell number.
- `snp` Required, matched_dir.
- `capture_virus` Required, matched_dir.

5th column:
- `dynaseq` Required, background snp file.

Example

Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.
```
$cat ./my.mapfile
fastq_prefix1	fastq_dir1	sample1
fastq_prefix2	fastq_dir2	sample1
fastq_prefix3	fastq_dir1	sample2

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
```

`--mod` Which type of script to generate, `sjm` or `shell`.

`--queue` Only works if the `--mod` selects `sjm`.

`--rm_files` Remove redundant fastq and bam files after running.

`--steps_run` Steps to run. Multiple Steps are separated by comma. For example, if you only want to run `barcode` and `cutadapt`, 
use `--steps_run barcode,cutadapt`.

`--outdir` Output directory.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--chemistry` Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:  
- `auto` Default value. Used for Singleron GEXSCOPE libraries >= scopeV2 and automatically detects the combinations.  
- `scopeV1` Used for legacy Singleron GEXSCOPE scopeV1 libraries.  
- `customized` Used for user defined combinations. You need to provide `pattern`, `whitelist` and `linker` at the 
same time.

`--pattern` The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  
- `C`: cell barcode  
- `L`: linker(common sequences)  
- `U`: UMI    
- `T`: poly T.

`--whitelist` Cell barcode whitelist file path, one cell barcode per line.

`--linker` Linker whitelist file path, one linker per line.

`--lowQual` Default 0. Bases in cell barcode and UMI whose phred value are lower than lowQual will be regarded as low-quality bases.

`--lowNum` The maximum allowed lowQual bases in cell barcode and UMI.

`--nopolyT` Outputs R1 reads without polyT.

`--noLinker` Outputs R1 reads without correct linker.

`--allowNoPolyT` Allow valid reads without polyT.

`--allowNoLinker` Allow valid reads without correct linker.

`--output_R1` Output valid R1 reads.

`--gzip` Output gzipped fastq files.

`--adapter_fasta` Addtional adapter fasta file.

`--minimum_length` Default `20`. Discard processed reads that are shorter than LENGTH.

`--nextseq_trim` Default `20`. Quality trimming of reads using two-color chemistry (NextSeq). 
Some Illumina instruments use a two-color chemistry to encode the four bases. 
This includes the NextSeq and the NovaSeq. 
In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
However, dark cycles also occur when sequencing “falls off” the end of the fragment.
The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.

`--overlap` Default `10`. Since Cutadapt allows partial matches between the read and the adapter sequence,
short matches can occur by chance, leading to erroneously trimmed bases. 
For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter. 
To reduce the number of falsely trimmed bases, the alignment algorithm requires that 
at least {overlap} bases match between adapter and read.

`--insert` Default `150`. Read2 insert length.

`--cutadapt_param` Other cutadapt parameters. For example, --cutadapt_param "-g AAA".

`--outFilterMatchNmin` Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.

`--out_unmapped` Output unmapped reads.

`--STAR_param` Other STAR parameters.

`--outFilterMultimapNmax` Default `1`. How many places are allowed to match a read at most.

`--starMem` Default `30`. Maximum memory that STAR can use.

`--gtf_type` Specify feature type in GTF annotation.

`--featureCounts_param` Other featureCounts parameters.

`--expected_cell_num` Default `3000`. Expected cell number.

`--cell_calling_method` Default `auto`. Cell calling methods. Choose from `auto` and `cellranger3`.

`--genomeDir` Required. Genome directory after running `celescope rna mkref`.

`--strand` gene strand file, the format is "geneID,+/-".

`--bg_cov` background snp depth filter, lower than bg_cov will be discarded. Only valid in csv format.

