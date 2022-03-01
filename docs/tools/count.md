## Features
- Cell-calling: Distinguish cell barcodes from background barcodes. 
- Generate expression matrix.
## Output
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
## Arguments
`--genomeDir` Required. Genome directory.

`--expected_cell_num` Default `3000`. Expected cell number.

`--cell_calling_method` Default `auto`. Cell calling methods. Choose from `auto` and `cellranger3`.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--bam` Required. BAM file from featureCounts.

`--force_cell_num` Default `None`. Force the cell number within (value * 0.9, value * 1.1).

