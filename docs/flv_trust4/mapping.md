## Features

- Extract candidate reads to assemble.

## Output
- `02.mapping/{sample}_bcrtcr.fq` All candidate reads(mapped to any V(D)J genes) sequence.
- `02.mapping/{sample}_bcrtcr_bc.fa` All candidate reads(mapped to any V(D)J genes) barcode.
- `02.mapping/{sample}_bcrtcr_umi.fa` All candidate reads(mapped to any V(D)J genes) umi.
- `02.mapping/{sample}_{chain}.fq` Candidate reads(mapped to {chain}) sequence. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.
- `02.mapping/{sample}_{chain}_bc.fa` Candidate reads(mapped to {chain}) barcode. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.
- `02.mapping/{sample}_{chain}_umi.fa` Candidate reads(mapped to {chain}) umi. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.
## Arguments
`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--match_fq2` R2 reads matched with scRNA-seq.

`--match_fq1` R1 reads matched with scRNA-seq.

`--ref` reference name.

`--seqtype` TCR/BCR seq data.

`--barcodeRange` Barcode range in fq1, INT INT CHAR.

`--umiRange` UMI range in fq1, INT INT CHAR.

