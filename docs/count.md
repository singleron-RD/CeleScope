# count

## Features
- Distinguish cells and background.
- Generate expression matrix.

## Input
- BAM file from step featureCounts.

## Output
- `{sample}_all_matrix` The expression matrix of all detected barcodes.

- `{sample}_matrix_10X` The expression matrix of the barcode that is judged to be the cell can be read in by calling the Seurat::Read10X function.

- `{sample}_matrix.tsv.gz` The expression matrix of the barcode that is judged to be the cell, separated by tabs.

- `{sample}_count_detail.txt.gz` There are four columns: barcode, gene ID, UMI count, read_count.

- `{sample}_counts.txt` Six columns:
	- Barcode: barcode sequence
	- readcount: reads num of each barcode
	- UMI2: UMI num (reads >= 2) for each barcode
	- UMI: UMI num for each barcode
	- geneID: Gene num for each barcode
	- mark UB: background CB: cell

- `{sample}_downsample.txt` Three columns：
	- percent: Percentage of sampled reads
	- median_geneNum
	- saturation

- `barcode_filter_magnitude.pdf` Barcode-UMI plot

## Parameters

`--bam` Required. BAM file from featureCounts.

`--force_cell_num` Default 'None'.

`--cells` Expect the number of cells, select the point with the largest slope within this value ±10% to determine the UMI threshold, see celescope.tools.count.barcode_filter_with_magnitud.

`--resue` Determine the cell based on the inflection point of the barcode-UMI curve, see celescope.tools.count.rescue_cells.

## Metrics
- Estimated Number of Cells : the number of barcodes considered as cell-associated.

- Fraction Reads in Cells : the fraction of uniquely-mapped-to-transcriptome reads with cell-associated barcodes.

- Mean Reads per Cell : the number of valid reads divided by the estimated number of cells.

- Median UMI per Cell : the median number of UMI counts per cell-associated barcode.

- Total Genes : the number of genes with at least one UMI count in any cell.

- Median Genes per Cell : the median number of genes detected per cell-associated barcode.

- Saturation : the fraction of UMI originating from an already-observed UMI.