# count

## Features
- Cell-calling: Distinguish cell barcodes from background barcodes. 
- Generate expression matrix.

## Input
- BAM file from step featureCounts.

## Output
- `{sample}_all_matrix` The expression matrix of all detected barcodes. Can be read in by calling the `Seurat::Read10X` function.

- `{sample}_matrix_10X` The expression matrix of the barcode that is identified to be the cell. Can be read in by calling the `Seurat::Read10X` function.

- `{sample}_matrix.tsv.gz` The expression matrix of the barcode that is identified to be the cell, separated by tabs.

- `{sample}_count_detail.txt.gz` 4 columns: 
	- barcode  
	- gene ID  
	- UMI count  
	- read_count  

- `{sample}_counts.txt` 6 columns:
	- Barcode: barcode sequence
	- readcount: reads num of each barcode
	- UMI2: UMI num (reads >= 2) for each barcode
	- UMI: UMI num for each barcode
	- geneID: Gene num for each barcode
	- mark UB: background CB: cell

- `{sample}_downsample.txt` 3 columns：
	- percent: percentage of sampled reads
	- median_geneNum: median gene number per cell
	- saturation: sequencing saturation

- `barcode_filter_magnitude.pdf` Barcode-UMI plot.

## Parameters

`--bam` Required. BAM file from featureCounts.

`--genomeDir` Required. Directory contains genome Fasta, GTF and refFLAT file. If this argument is not provided, you need to provide `--gtf`.

`--gtf` GTF file path.

`--expected_cell_num` Default `3000`. Expected cell number.

`--cell_calling_method` Default `auto`. Cell calling methods.Choose from `auto`, `cellranger3` and `inflection`.

`--force_cell_num` Default `None`. Force the cell number to be this value ± 10%.

## Metrics
- Estimated Number of Cells : the number of barcodes considered as cell-associated.

- Fraction Reads in Cells : the fraction of uniquely-mapped-to-transcriptome reads with cell-associated barcodes.

- Mean Reads per Cell : the number of valid reads divided by the estimated number of cells.

- Median UMI per Cell : the median number of UMI counts per cell-associated barcode.

- Total Genes : the number of genes with at least one UMI count in any cell.

- Median Genes per Cell : the median number of genes detected per cell-associated barcode.

- Saturation : the fraction of UMI originating from an already-observed UMI.