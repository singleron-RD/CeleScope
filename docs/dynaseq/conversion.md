## Features
- Get conversion pos in each read.
    - Get snp info. 

## Output
- `{sample}.PosTag.bam` Bam file with conversion info.
- `{sample}.PosTag.csv` SNP info in csv format.
## Arguments
`--strand` gene strand file, the format is "geneID,+/-".

`--bam` featureCount bam(sortedByCoord), must have "MD" tag, set in star step.

`--cell` barcode cell list.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

