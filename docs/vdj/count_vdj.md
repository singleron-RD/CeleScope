## Features
- Cell-calling based on barcode-UMI rank.    
- Summarize clonetypes infomation.

## Output
- `{sample}_cell_confident.tsv` The clone type of VDJ cell barcode, each chain occupies one line.

- `{sample}_cell_confident_count.tsv` The clone type of VDJ cell barcode, each cell occupies one line.

- `{sample}_clonetypes.tsv` The count and percentage of each clonetypes of VDJ cell barcode.

- `{sample}_match_clonetypes.tsv` When summarize clonetypes, only consider barcodes in the match scRNA-Seq library. 
This file will only be produced when the `match_dir` parameter is provided.
## Arguments
`--type` Required. `TCR` or `BCR`.

`--UMI_min` Default `auto`. Minimum UMI number to filter. The barcode with UMI>=UMI_min is considered to be cell.

`--iUMI` Default `1`. Minimum number of UMI of identical receptor type and CDR3. 
For each (barcode, chain) combination, only UMI>=iUMI is considered valid.

`--UMI_count_filter_file` Required. File from step mapping_vdj.

`--match_dir` Match celescope scRNA-Seq directory.

`--matrix_dir` Match celescope scRNA-Seq matrix directory.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

