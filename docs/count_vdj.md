# count_vdj

## Features
- Cell calling.
- Calculate clonetypes.

## Input
- `{sample}_UMI_count_filtered1.tsv` from mapping_vdj.

## Output
- `{sample}_cell_confident.tsv` The clone type of VDJ cell barcode, each chain occupies one line.

- `{sample}_cell_confident_count.tsv` The clone type of VDJ cell barcode, each cell occupies one line.

- `{sample}_clonetypes.tsv` The count and percentage of each clonotype of VDJ cell barcode.

- `{sample}_match_clonetypes.tsv` The count and percentage of each clonotype at the intersection of VDJ cell barcode and sc-RNA-Seq cell barcode. The file will only be produced when the match_dir parameter is provided.

## Parameters

`--type` Required. `TCR` or `BCR`.

`--UMI_count_filter1_file` Required. File from step mapping_vdj.

`--match_dir` Default `None`. CeleScope rna analysis path. 

`--UMI_min` Default `auto`. Minimum UMI number to filter. The barcode with UMI>=UMI_min is considered to be cell.

`--iUMI` Default `1`. Minimum number of UMI of identical receptor type and CDR3. For each (barcode, chain) combination, only UMI>=iUMI is considered confident. 

## Metrics
- Estimated Number of cells : number of barcodes considered as cell-associated.

- Median UMIs per Cell : median number of UMI mapped to each chain per cell.

- Cell with Heavy and Light Chain : cells with as least 5 UMI mapped to each chain.

- Cell with Barcode Match : cells with barcode matched with scRNA-seq library.

- Cell with barcode Match, Heavy and Light chain : cell with matched barcode and with as least 5 UMI mapped to each chain.


