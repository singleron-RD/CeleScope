## Features
- Align R2 reads to IGMT(http://www.imgt.org/) database sequences with blast.

## Output
- `{sample}_airr.tsv` The alignment result of each UMI.
A tab-delimited file compliant with the AIRR Rearrangement schema(https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

- `{sample}_UMI_count_unfiltered.tsv` UMI reading for each (barcode, chain, VJ_pair) combination.

- `{sample}_UMI_count_filtered.tsv` For each (barcode, chain) combination, only the record with the 
most VJ_pair UMI reads is kept.
## Arguments
`--ref_path` reference path for igblast.

`--type` TCR or BCR.

`--species` Default human. human or mouse.

`--fasta` Required. Input fasta file.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

