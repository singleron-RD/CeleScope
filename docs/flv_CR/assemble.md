## Features

- TCR/BCR Assemble by Cellranger.

- Generate Mapping, Cells, V(D)J annotations metrics in html.

## Output

- `03.assemble/{sample}` Cellranger vdj results.
## Arguments
`--ref_path` reference path for cellranger.

`--soft_path` soft path for cellranger.

`--other_param` Other cellranger parameters.

`--mem` memory(G).

`--seqtype` TCR or BCR.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fqs_dir` fastq dir after convert.

