# Change Log
 
## [Unreleased] - 2021-05-17

 
### Added

- Assay `rna` outputs .h5 file in 06.analysis directory.

### Changed

- Update Seurat from 2.3.4 to 4.0.1.

### Fixed

### Removed

- Assay `rna` no longer outputs tab-delimited expression matrix file in 05.count directory.

 
## [1.1.9] - 2021-04-25
 
### Added

- Add parameter `--coefficient`  to `celescope tag count_tag` and `multi_tag`
    
    Default `0.1`. Minimum signal-to-noise ratio is calulated as `SNR_min = max(median(SNRs) * coefficient, 2)`

- Add `.metrics.json`
 
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
