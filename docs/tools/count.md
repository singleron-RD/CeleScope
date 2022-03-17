## Features
- Cell-calling: Distinguish cell barcodes from background barcodes. 
- Generate expression matrix.
## Output
- `{sample}_raw_feature_bc_matrix` The expression matrix of all detected barcodes in [Matrix Market Exchange Formats](
    https://math.nist.gov/MatrixMarket/formats.html). 
- `{sample}_filtered_feature_bc_matrix` The expression matrix of cell barcodes in Matrix Market Exchange Formats. 
- `{sample}_count_detail.txt.gz` 4 columns: 
    - barcode  
    - gene ID  
    - UMI count  
    - read_count  
- `{sample}_counts.txt` 6 columns:
    - Barcode: barcode sequence
    - readcount: read count of each barcode
    - UMI2: read count with reads per UMI >= 2 for each barcode
    - UMI: UMI count for each barcode
    - geneID: gene count for each barcode
    - mark: cell barcode or backgound barcode.
        `CB` cell  
        `UB` background  
- `{sample}_downsample.tsv` Subset a fraction of reads and calculate median gene number and sequencing saturation.
## Arguments
`--genomeDir` Required. Genome directory.

`--expected_cell_num` Default `3000`. Expected cell number.

`--cell_calling_method` Default `auto`. Cell calling methods. Choose from `auto` and `cellranger3`.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--bam` Required. BAM file from featureCounts.

`--force_cell_num` Default `None`. Force the cell number to be this number.

