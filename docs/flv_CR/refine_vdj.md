## Features

- Refine barcodes where "is_cell=False" and have multi productive chains.

- There are three methods to filter noise: SNR, AUTO, NOT_FILTER
## Arguments
`--not_refine` Run the VDJ refine step.

`--filter_method` filter noise method.

`--coeff` coefficient will affect auto and snr noise filter, recommend 1.5 for auto, 10 for snr.

`--seqtype` TCR or BCR.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--assemble_out` directory of cellranger assemble result.

`--match_dir` scRNA-seq match directory.

`--barcode_convert_json` json file.

