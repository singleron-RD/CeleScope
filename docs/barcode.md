# barcode

## Major function
- Barcode deumltiplexing and filtering.

## Input
- Paired-end FASTQ files.

## Output
- `{sample}_2.fq` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of the read name is `{barcode}_{UMI}_{read ID}`.

- `fastqc.zip` `fastqc.html` Fastqc results.

## Parameters

`--fq1` Required, FASTQ R1 reads. Multiple FASTQ files are seperated by comma. 

`--fq2` Required, FASTQ R2 reads. Multiple FASTQ files are seperated by comma. 

`--chemistry` Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:  
`auto` Default value. For Singleron GEXSCOPE libraries >= scopeV2, automatically detects the used chemistry.  
`scopeV1` For legacy Singleron GEXSCOPE scopeV1 libraries.  
`customized` For user defined combinations. You need to provide `pattern`, `whitelist` and `linker` at the same time.

`--pattern` the pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number of bases.  
`C`: cell barcode  
`L`: linker(common sequences)  
`U`: UMI    
`T`: poly T

`--whitelist` cell barcode whitelist file path, one cell barcode per line.

`--linker` linker whitelist file path, one linker per line.

`--lowQual` Bases in cell barcode and UMI whose phred value are lower than lowQual will be regarded as low-quality bases. Default=0.  

`--lowNum` The maximum allowed lowQual bases in cell barcode and UMI. Default=2.

`--nopolyT` Outputs R1 reads without polyT.
    
`--noLinker` Outputs R1 reads without correct linker.

`--allowNoPolyT` Allow reads without polyT.

`--allowNoLinker` Allow reads with incorrect linker.

## Metrics
-  Raw Reads: total read pairs from FASTQ files.

-  Valid Reads: reads pass filtering(filtered: reads without poly T, reads without linker, reads without correct barcode or low quality reads).

-  Q30 of Barcodes: percent of barcode base pairs with quality scores over Q30.

-  Q30 of UMIs: percent of UMI base pairs with quality scores over Q30.