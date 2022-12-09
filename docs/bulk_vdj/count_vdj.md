## Features  
- Summarize clonetypes infomation.

## Output
- `{sample}_filtered_annotations.csv` Annotations for each CDR3.

- `{sample}_clonetypes.csv` The count and percentage of each CDR3 of reads.
## Arguments
`--type` Required. `TCR` or `BCR`.

`--productive_file` Required. File from step mapping_vdj.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

