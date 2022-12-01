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
 - If two or more groups of reads have the same barcode and UMI, but different gene annotations, the gene annotation with the most supporting reads is kept for UMI counting and the other read groups are discarded. In case of a tie for maximal read support, all read groups are discarded. In previous versions, both read groups were preserved.
 
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
 - Add a sub-command `celescope utils mkgtf` which is similar to `cellranger mkgtf`. After using this command, only the lines with gene_biotype as protein_coding, lncRNA, antisense and VDJ related genes will be kept in gtf. This removes 2 mitochondrial ribosomal RNAs (mt-rRNA) and 22 mitochondrial transfer RNAs (mt-tRNA) from gtf. The detected mitochondrial gene UMI is decreased.

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
 - Rename the cell-calling method from `cellranger3` to `EmptyDrops_CR`. Make `EmptyDrops_CR` the default method.
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

- Add parameter `--cell_calling_method` to `celescope rna count` and `multi_rna`.

    `--cell_calling_method` can be one of:  
    - `auto`: Same result as v1.1.7.  
    - `cellranger3`: Refer to the cell_calling algorithm of cellranger3, and the result is similar to cellranger3.  
    - `reflection`: Use the inflection point of the barcode-rank curve as the UMI threshold. The minimum UMI value is changed from initial threshold / 10 to initial threshold / 2 to prevent the use of a lower inflection point when there are multiple inflection points.  

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




