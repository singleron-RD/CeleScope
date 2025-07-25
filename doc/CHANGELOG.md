## [2.7.4] - 2025-07-14
- Fixed: In `bulk_rna split_fastq`, the sequences are now the original sequences from the FASTQ file. Previously, if a read was mapped to the reverse strand, the split FASTQ contained the reverse-complemented sequence.

## [2.7.3] - 2025-07-09
- Changed: The method used for automatic chemistry detection (`--chemistry auto`) has been updated. Previously, the 6 bp linker sequence following the UMI was used to distinguish between GEXSCOPE-V1 and flv_rna. This version now uses the presence of a 'C' at the 57th base instead, providing more robust classification in cases where the 6 bp linker may be affected by sequencing errors.

## [2.7.2] - 2025-07-09
 - Changed: In `bulk_vdj` analysis, retain only UMIs whose CDR3 amino acid sequence is longer than 5 residues and starts with 'C'.

## [2.7.1] - 2025-07-08
 - Added: `--outSAMtype None` for no BAM output.

## [2.7.0] - 2025-07-04
 - Added: support for `bulk_rna-V3` auto detection.
 - Added: well metrics in `bulk_vdj` HTML report.
 - Changed: introduced the required `--well_sample` parameter in `bulk_vdj`. 
 - Changed: optimized the consensus algorithm so that leading Ns no longer affect the downstream base calls. This improvement can increase the V(D)J mapping rate after consensus.
 - Removed: `cutadapt` step in `vdj` and `bulk_vdj`.

## [2.6.1] - 2025-06-10
 - Added: log output support for shell mode.

## [2.6.0] - 2025-05-29
 - Added: `--report_soloFeature` in `starsolo`. If you do not want to include intron reads in the analysis, use `--report_soloFeature Gene`.
 - Added: parameters in the HTML report.

## [2.5.0] - 2025-05-13
 - Added: BAM splitting functionality and read-level metrics in `bulk_rna`.
 - Added: `--limitBAMsortRAM` parameter in `starsolo`.

## [2.4.0] - 2025-04-02
 - Added: support for `flv_rna-V2` chemistry.
 - Fixed: "customized" option in chemistry was not working.

## [2.3.2] - 2025-03-11
 - Fixed: `bulk_rna-V2` was auto-detected as `bulk_rna-V1`.

## [2.3.1] - 2025-03-06
 - Changed: Improve the split_fastq speed of `bulk_rna`.

## [2.3.0] - 2025-03-04
 - Changed: `bulk_rna` added automatic chemistry detection and FASTQ splitting functionality, replaced featureCounts with STARsolo for quantification, introduced the required `--well_sample` parameter.

## [2.2.2] - 2025-02-12
 - Fixed: Pin plotly version == 5.24.1 to avoid compatibility issues.
 
## [2.2.0] - 2025-02-05
 - Added: pathSeq workflow for single-cell 16S.
 - Added: support for species other than human or mouse in `flv_trust4`.

## [2.1.0] - 2024-09-10
 - Added: support for `rna_5p3p`.
 - Fixed: Pin numpy version==1.26.0.
 - Fixed: Wrong delimiter results in empty filter_contig.fasta.
 - Changed: When some read lengths in R1 read are shorter than required, no error is reported.

## [2.0.7] - 2023-12-01
 ### `rna` and `dynaseq`
  - Adjust the denominator of `Reads Assigned To {Exonic, Intronic, Intergenic, Antisense} regions` from valid reads to mapped reads. This way Exonic + Intronic + Intergenic + Antisense = 100%

 ### General improvments
 - Bug fix: Introns was added to gtf when running `celescope utils mkgtf`. Otherwise the featureCounts step will report an error for not finding intron regions.
 - Removed the --gzip argument.

## [2.0.6] - 2023-11-16
 ### General improvments
 - Add support for legacy chemistry `scopeV1`.

## [2.0.5] - 2023-11-15
 - Add a new panel `CHIP`.
 - Merge all overlapping regions in bed files. Fix KRAS gene in lung_1.bed.

## [2.0.4] - 2023-11-07
 ### General improvments
 - Fix a bug in parsing matrix files from `matched_dir`.

## [2.0.3] - 2023-10-31
 ### `tag`, `citeseq` and `sweetseq`
 - Remove cutadapt from steps.

 ### `citeseq`
 - Remove low expression antibody from HTML Umap plot.

## [2.0.2] - 2023-10-30
 ### `rna` and `dynaseq`
 - `celescope rna cells` add a new argument `--soloCellFilter` for cell filtering of the raw count matrix.

## [2.0.1] - 2023-10-27
 ### `rna` and `dynaseq`
 - Add a metric in HTML report - Corrected Barcodes.
 - `celescope rna cells` add a new argument `--min_gene` to filter cells with low number of genes.

## [2.0.0] - 2023-10-25
 ### `rna` and `dynaseq`
 - Use Starsolo to quantify the gene expression.

 - The `celescope rna cells` step has been added, which can perform mitochondrial quality control and forced cell number on the existing results, and generate the corresponding expression matrix and HTML report. The original expression matrix and HTML report will be prefixed with 'default_'; the cells_summary in .metrics.json and .data.json will also be prefixed with 'default_'.

 - Update the STAR version from 2.6.1 to 2.7.11a; since the STAR index file is incompatible between the two versions, the genome needs to re-run mkref. scanpy upgraded to 1.9.3.

 - When using celescope utils mkgtf to filter gtf, if the program finds that there is an MT chromosome in gtf, an additional mt_gene_list.txt file will be generated, including the gene_name of the gene on the MT chromosome. This file can be provided to the mkref command during mkref, so that the analysis step can count the percentages of these genes. This file will also be used in subsequent mitochondrial quality control (celescope rna cells).

 - The main output files (bam, expression matrix, h5ad, etc.) are now placed in the `{sample}/outs/` directory. Expression matrix files use gzip compression.
 
 - Added HTML reporting metrics

  `Q30 of RNA Reads` Fraction of RNA read bases with quality score >= 30.

  `Reads Mapped Uniquely To Transcriptome` Reads that mapped to a unique gene in the transcriptome. These reads are used for UMI counting.

  `Reads Assigned Antisense To Gene` Reads that assigned to the opposite strand of genes.

  `Mean Used Reads per Cell` The number of uniquely-mapped-to-transcriptome reads per cell-associated barcode.

 - In bam and matrix, the three barcode segments are now separated by underscore '_'.

## [1.17.0] - 2023-08-11
 ### `bulk_vdj`
 - Add support for AccuraCode VDJ(bulk vdj).
 
 ### `dynaseq`
 - Various updates. https://github.com/singleron-RD/CeleScope/pull/252
 
 ### `tag`
 - Fix a bug in the `split_tag` step.


## [1.16.1] - 2023-07-07
 ### `rna` and `dynaseq`
 - Reduce the memory consumption of `count`; improve speed of `featureCounts`.

## [1.16.0] - 2023-07-05
 ### `rna` and `dynaseq`
 - Allow all numeric "gene_id" in gtf files. In previous versions, all numeric "gene_id" would cause an error in the `count` step.
 - Merge `barcode`, `cutadapt`, `STAR` into one step `prep_map`. Intermediate fastqs are no longer output.
 - Change the algorithm of sequencing saturation: Calculate the percentage of PCR duplicate reads instead of the percentage of duplicate UMI. 

 ### `snp`
 - When annotating human variants, [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) is used by default instead of the original ensembl 99. This helps to improve annotation quality.

 ### `tag`
- When using `--split_bam`, reads that map to intergenic regions will also be written to each split bam.

 ### General improvements
 - fix #127.
 - Remove redundant conda and python packages to improve installation speed.
 - Enforce a maximum number of threads of 20. The number of threads greater than 20 will not bring benefits during STAR mapping, and may cause disk errors during samtools sort bam.

## [1.15.2] - 2023-05-24
 ### `flv_CR`
 - Fixed a bug when running the refine step on BCR data.

## [1.15.1] - 2023-05-22
 ### `flv_CR`
 - An additional refine step is added to rescue some cells from cells assembling multiple chains (possibly due to background contamination).

 ### `vdj`
 - Allows a chain to be missing in all cells (e.g. all cells have no IGL chain); previous versions reported an error.
 - When reporting UMI number mapping to each chain, use mean instead of the previous median. This makes it easier to understand in the presence of multiple light chains.

 ### General improvments
 - Reduce the size of HTML reports.

## [1.15.0] - 2023-04-11
 ### `snp`
 - Added image output of Otsu's method and genotype mapping on UMAP graph.
 - Updated genes in the panel.

 ### `citeseq`
 - Fixed the problem that the T-SNE graph was shrunk when the protein name was too long.

 ### `flv_CR` and `flv_trust4`
 - Fixed the problem that the color of barcode rank plot in HTML report could not be displayed correctly when the cell density is low.

 ### General improvments

 - Fixed the version of pandas to 1.4.2 to avoid compatibility issues caused by pandas 2.0.0.
 - Change the default value of STAR --outFilterMatchNmin from 0 to 50 to avoid the interference of short fragments (especially probe sequences in FocuScope) on gene quantification.
 - Multiple barcode whitelists are allowed.

## [1.14.1] - 2023-01-11
 ### `snp`
 - Add an argument `--database`. Data from non-human species can also be annotated with this argument.
 - Change the default filtering method from `auto` to `otsu` to improve sensitivity.
 - Fix an issue that could result in blank gene annotations.
 - Fix an issue where amino acid names were not displaying correctly.

## [1.14.0] - 2022-12-20
 ### `rna`
 - Revert changes to the rna pipeline made in version 1.13.0. 
 
 > If two or more groups of reads have the same barcode and UMI, but different gene annotations, the gene annotation with the most supporting reads is kept for UMI counting and the other read groups are discarded. In case of a tie for maximal read support, all read groups are discarded.

 ### `snp`
 - Replace ANNOVAR with SnpEff.

 ### `citeseq`
 - Optimize the display of images in HTML reports.

## [1.14.0b0] - 2022-12-01
 ### `vdj`
 - Replace Mixcr with IgBlast.

 ### `snp`
 - When using `otsu` or `auto` , two additional filtering conditions are added: Minimum variant allele frequency and Minimum supporting reads . This will remove some low-confidence heterozygous variants that cannot be filtered out by the `otsu` or `auto` algorithm.
 - Some Splits reads (e.g. spanning splicing events in RNAseq data) will result in a small number of variants that are not in the target gene list. These variants are removed.

## [1.13.0] - 2022-10-28

 ### `sweetseq`
 - Added support for single cell data generated with ProMoSCOPE<sup>TM</sup> Single Cell Glycosylation Detection Kits.

 ### `capture_virus`
 - The `capture_virus` pipeline now supports adding gtf as an input file to generate an expression matrix of viral genes.

 ### `rna`
 - If two or more groups of reads have the same barcode and UMI, but different gene annotations, the gene annotation with the most supporting reads is kept for UMI counting and the other read groups are discarded. In case of a tie for maximal read support, all read groups are discarded. In previous versions, both read groups were kept.
 
 ### General improvments
 - There are some gtf missing lines with annotation as gene. Such gtf will report an error when running featureCounts. The `celescope utils mkgtf` command can now add missing lines for such gtf.

## [1.12.1] - 2022-09-16
 ### General improvments
 - Fixed version of package `matplotlib`. The latest version of matplotlib 3.6.0 causes the error: `ModuleNotFoundError: No module named 'matplotlib._contour'`.

## [1.12.0] - 2022-09-08
 ### `rna`
 - Introns will be included in the `rna` analysis by default.

 ### `vdj`
 - Add an argument `--mixcr_mem` to avoid the `Invalid maximum heap size` error reported by mixcr. https://github.com/milaboratory/mixcr/issues/588

 ### General improvments
 - `Reads without poly T` are not filtered by default. Remove the argument `--allowNoPolyT` and add a new argument `--filterNoPolyT`.
 - Add a sub-command `celescope utils mkgtf`. After using this command, only the lines with gene_biotype as protein_coding, lncRNA, antisense and VDJ related genes will be kept in gtf. This removes 2 mitochondrial ribosomal RNAs (mt-rRNA) and 22 mitochondrial transfer RNAs (mt-tRNA) from gtf. The detected mitochondrial gene UMI is decreased.

## [1.11.1] - 2022-08-10
 ### General improvments
 - When making a STAR index, `--genomeSAindexNbases` can be automatically inferred based on genome size.
 - Fix an issue that `Reads Assigned To Intronic Regions` could have negative values.
 - Fix the problem that the excel download button of the marker gene list was removed.

 ### `snp`
 - Fix an issue where Annovar could generate extra lines in multianno.txt file.

## [1.11.0] - 2022-07-07
 ### `rna`
 - All genes from the gtf file will be written to `genes.tsv`. https://github.com/singleron-RD/CeleScope/issues/81
 - The `celescope rna mkref` command no longer generates refflat files; removes the number of bases mapped to each region; adds the number of reads mapped to each region.
 
 ### `flv_CR`
 - fix a bug that `VJ Spanning (IGH, IGL) pair` will cause an KeyError if the order of [IGH, IGL] is not satisfied. https://github.com/singleron-RD/CeleScope/pull/145 


## [1.11.0b0] - 2022-06-21
 ### `flv_vdj`
 - Chemistry `flv_rna`(full length VDJ matched mRNA) and `flv` can be auto-detected. 
 - Remove assay `flv_trust4_split`.

 ### `snp`,`capture_virus` and `fusion`
 - Change the default `--umi_threshold_method` from `auto` to `otsu`.

 ### `vdj`
 - Use an improved cell-calling method that only match cell barcodes are considered. 
 - Change how median IGH, IGK and IGL UMIs are calculated. Now all cells(include 0 UMI count cells) are taken into consideration.
 - Add cell with chain pair metrics, e.g. "Cells with (IGH, IGK) pair".

 ### `mkref`
 - Change parameters for generating refFlat files to tolerate non-ensembl-compliant gtf.

 

## [1.10.0] - 2022-04-22
 ### `fl_vdj`
 - Add 3 new assays: `flv_CR`, `flv_trust4` and `flv_trust4_split`.

 ### `rna` and `dynaseq`
 - Limit the marker genes of each cluster in the HTML report to a maximum of 50 to avoid the slow opening of the report.
 - Add read_saturation to the HTML report and `{sample}/*.count/downsample.txt`

 ### `snp`
 - Add panel `blood_1`.

## [1.9.0] - 2022-04-01
 ### `rna` and `dynaseq`
 - Fix an issue that mitochondrial percent is not added to metrics.

 ### `snp`,`capture_virus` and `fusion`
 - When calculating the `auto` threshold, the default coefficient changes from 10 to 3. This will make the filtering more stringent.

## [1.8.1] - 2022-03-23
 ### General improvments
 - Fix an issue where the matrix suffix `filtered_feature_bc_matrix` introduced in v1.8.0 is not recognized when parsing match_dir.

## [1.8.0] - 2022-03-17
 ### `rna` and `dynaseq`
 - Replace `Seurat` with `scanpy`.
 - Add read_saturation to downsample file.

 ### `snp`,`capture_virus` and `fusion`
 - When calculating `otsu` threshold, use `math.ceil` instead of `int`.
 
 ### General improvments
 - Fix an issue where the conditions for detecting scopeV2.0.1 are too loose. (#108)
 - Move `sjm.job` from `./log/` to `./sjm/`.
 - Change output file suffix.
    - raw_matrix: `all_matrix` -> `raw_feature_bc_matrix`
    - fitered_matrix: `matrix_10X`-> `filtered_feature_bc_matrix`


## [1.7.2] - 2022-02-10

 ### `vdj`
 - "Cell with Barcode Match, TRA and TRB": When calculating the percentage, the denominator is `Cell with Barcode Match`. The denominator used previously was `Estimated Mumber of Cells`.

 ### `capture_virus` and `fusion`
 - Make UMI correction optional. If you do not want to perform UMI correction, use `--not_correct_UMI`.
 - Add a filtering step to filter UMI with supporting reads less then read_threshold. 

 ### `dyanseq`
 - modify vcf base [#96](https://github.com/singleron-RD/CeleScope/pull/96).

 ### General improvments
 - Remove the redundant `--assay` parameter.
 - Add the `--queue` argument for `sjm` job submission.

## [1.7.1] - 2022-01-17

 ### `rna`
 - Fix a bug with mt_gene_list (#92)

 ### `snp`
  - Add a fitering step: `celescope snp filter_snp` with arguments `--threshold_method`. Choices can be one of 
    - `auto` : Default method. Using a method similar to cell calling method.
    - `otsu` : Counts are first log transformed and then the threshold is determined by [Otsu's method](https://en.wikipedia.org/wiki/Otsu%27s_method).
    - `hard` : Using user provided UMI threshold.
    - `none` : Do not perform filtering. 

## [1.7.0] - 2021-12-28
 ### `capture_virus`
  - Add documents.
  - Add UMI filtering options - 'auto'.
  - Add UMI correction.
  - Extract the filtering process into a separate step.
  - Remove the `otsu` choice from `--min_support_read`.

 ### `tag`
  - Remove `--marker_file` argument from `analysis_tag`.
 
  ### General improvments
  - Optimize the display of t-SNE plots in HTML reports. 

## [1.6.1] - 2021-12-01
 ### `snp`
  - Increase `max-depth` from 1M to 100M to avoid calculation errors at very deep sites.

 ### General improvments
  - Fix a bug that `featureCounts` didn't output name sorted bam file.
  - Fix a bug that `Median Enriched Reads per Valid Cell` in `target_metrics` used all cell number as denominator, not valid cell number.
  - Fix a problem that the Q30 metrics in the html report showed too much precision.

## [1.6.0] - 2021-11-29

 ### `snp` 
   - Improve the speed and memory usage.
   - Format the output files to avoid redundant information.
   - Fix an issue that the number of cells of each genotypes were not calculated correctly.
   - Remove `CID` and `VID` in all the output files.   

 ### General improvments
  - Improve the speed of `add_tag`.
  - Remove support for deprecated cell-calling method `inflection`.

## [1.5.2] - 2021-11-02
 ### General improvments
 - Add auto-detection for chemistry `scopeV3.0.1`.
 - Remove support for deprecated chemistry `scopeV2.0.0` and `scopeV2.1.0`.

## [1.5.1] - 2021-10-28
 ### `snp`
 - Add `--panel` option.
 - Add output files.
 - Improve the speed of `variant_calling`.


 ### `fusion`
 - Fix a bug where `celescope.fusion.count_fusion` do not recognize fusion reads correctly.

 ### `dynascope`
 - Fix some bugs (#62).

 ### General improvments
  - Use local css instead of urls to avoid connection issues.


## [1.5.0] - 2021-09-09
 ### Added

 - Add `--split_vdj` in `celescope.tag.split_tag` to split vdj library according to tag assignment.
 - Add Docs for `snp`.

 ### Changed

 - `snp.variant_calling` change `--do-not-fix-overhangs` from `false` to `true` to avoid omitting variants in the overhang region.

 ### Fixed

 - Fix a memory leak in `snp.variant_calling`.

 ### Removed

## [1.4.0] - 2021-08-24

 ### Added

 - Add `min_consensus_read` argument to `celescope.tools.consensus`. If 
    1. the percentage of the most common residue type > threshold;
    2. most common residue reads >= min_consensus_read;
    
then we will add that residue type, otherwise an ambiguous character will be added.

 ### Changed
 - By default, use otsu method to calculate `min_support_read` for `capture_virus`.

 - By default, use otsu method to calculate `min_support_read` for `snp`.

 - Improved the code of `celescope.snp.variant_calling` to reduce memory usage and speed up analysis.

 - Add serveral arguments in `vdj` and `tag` to support WDL workflow.

 ### Fixed

 - Fix a bug in `celescope.tag.count_tag` when there is no `Undetermined` cells.

 - Fix a bug in `celescope.snp.variant_calling` when there are multiple variants at the same loci in a cell.

 ### Removed

## [1.3.2] - 2021-07-09
### Added

- Add `dynascope` assay.
- Add `celescope tag split_tag`. 

### Changed
- Change fastq file pattern of mapfile: Remove * before library_id.
- `celescope.tools.count_capture_virus`: Change `min_support_read` from 1 to 2.

### Fixed
- `celescope.tools.count` will report an error when there are multiple gtf or refFlat file under `genomeDir`.

### Removed
- `celescope.tools.utils.glob_genomeDir`

## [1.3.1] - 2021-06-09
### Added

- Add wdl workflow.

- Add Seurat hashtag method in `celescope tag count_tag`. To get Seurat hashtag output, use `--debug`. However, there was a unsolved problem with this method: https://github.com/satijalab/seurat/issues/2549.

### Changed

- `{sample}_UMI_count_filtered1.tsv` in mapping_vdj changed to `{sample}_UMI_count_filtered.tsv` (remove `1` after filtered)

### Fixed and Removed

- Remove h5 file generation in R to avoid memory issues.


## [1.3.0] - 2021-05-28
 
### Added

- `mkref` subcommand. See `celescope rna mkref`, `celescope fusion mkref` and `celescope virus mkref` for details.

### Changed

- Change the way to handle duplicate gene_name and gene_id in gtf file.

Previous:

    - one gene_name with multiple gene_id: "_{count}" will be added to gene_name.
    - one gene_id with multiple gene_name: newer gene_name will overwrite older gene_name.
    - duplicated (gene_name, gene_id): "_{count}" will be added to gene_name.

Now:

    - one gene_name with multiple gene_id: "_{count}" will be added to gene_name.
    - one gene_id with multiple gene_name: error.
    - duplicated (gene_name, gene_id): ignore duplicated records and print a warning.

### Fixed

- Fix `count tag` metrics order in merge.xls

### Removed

- Remove `--fusion_pos` from `celescope.fusion.count_fusion`

 
## [1.2.0] - 2021-05-19
 
### Added

- Assay `rna` outputs .h5 file in 06.analysis directory.

### Changed

- Update Seurat from 2.3.4 to 4.0.1.

- `--genomeDir` in `celescope.fusion.star_fusion` changed to `--fusion_genomeDir` to avoid misunderstanding.

- Step `star` sort bam by samtools instead of STAR to avoid potential `not enough memory for BAM sorting` error: https://github.com/alexdobin/STAR/issues/1136

### Removed

- Assay `rna` no longer outputs tab-delimited expression matrix file in 05.count directory.

 
## [1.1.9] - 2021-04-25
 
### Added

- Add parameter `--coefficient`  to `celescope tag count_tag` and `multi_tag`
    
    Default `0.1`. Minimum signal-to-noise ratio is calulated as `SNR_min = max(median(SNRs) * coefficient, 2)`

- Add `.metrics.json`

- Add `scopeV1` chemistry support.
 
### Changed
  
- Optimize speed and memory usage of step `barcode`(~2X faster) and `celescope.tools.count.downsample`(~15-25X faster, 1/2 memory usage).

- Change filtering of linker from allowing two mismatches in total to two mismatches per segment; this will slightly increase the valid reads percentage.

- Default output fastq files of `barcode` and `cutadapt` are not gzipped. Use `--gzipped` to get gzipped output.

- Change the display of Barcode-rank plot in html report. 

### Fixed

- Fix a bug that `celescope.tools.barcode.mismatch` cannot output all sequences correctly when n_mismatch>=2.

- Fix an error when Numpy >= 1.2.0.

- VDJ merge.xls can display all the metrics correctly.

### Removed

- Remove fastqc from `barcode` step.

 
## [1.1.8] - 2021-03-26
 
### Added

- Add read consensus to VDJ pipeline. 

    A consensus step was added before mapping to merge all the reads of the same
    (barcode, UMI) into one UMI. For defailed consensus algorithm, refer to `celescope.tools.consensus`.  
    multi_vdj adds the parameter `--not_consensus` that you can skip the consensus step, and get the same results as v1.1.7.   

- Add parameter `--species` to `celescope vdj mapping_vdj` and `multi_vdj`.

    `--species` can be one of:
    - `hs`: human
    - `mmu`: mouse

- Add 4 tags to featureCounts bam.

    - `CB`: cell barcode
    - `UB`: UMI
    - `GN`: gene name
    - `GX`: gene id

- Add `--STAR_param` to `celescope rna STAR`

    Additional parameters of STAR can be passed into the `STAR` step.

### Changed

- One sample can have different chemistry fastq in mapfile.  Version <= v1.1.7 will report this as an error.

- Gtf file can be gzipped.

- `multi_rna` can use 3 paramters: `--STAR_index`, `--gtf` and `--refFlat` instead of `--genomeDir` 

- Step `snpCalling` use mutract.


## [1.1.7] - 2020-12-16

### Added

- Automatically detect Singleron chemistry version.

### Changed

- FeatureCounts use strand specificity.

- Cutadapt default `overlap` change from `5` to `10`.

- VDJ sort `NA` last.

- `match clonetypes` are sorted by barcode_count(Frequency) first, then clonetype_ID.




