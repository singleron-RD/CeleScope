# mapping_tag

## Features
- Align R2 reads to the tag barcode fasta.

## Input
- R2 reads from step cutadapt.

## Output

- `{sample}_read_count.tsv` tab-delimited text file with 4 columns.

    `barcode` cell barcode  
    `tag_name`  tag name in barcode_fasta  
    `UMI`   UMI sequence  
    `read_count` read count per UMI  

## Parameters

`--fq_pattern` Required. R2 read pattern. The number after the letter represents the number of bases. 

`L` linker(common sequences)  
`C` tag barcode  


`--barcode_fasta` Required. Tag barcode fasta file. It will check the mismatches between tag barcode sequence in R2 reads with all tag barcode sequence in barcode_fasta. It will assign read to the tag with mismatch < len(tag barcode) / 10 + 1. If no such tag exists, the read is classified as invalid.
```
>tag_0
GGGCGTCTGTGACCGCGTGATACTGCATTGTAGACCGCCCAACTC
>tag_1
TTCCTCCAGAGGAGACCGAGCCGGTCAATTCAGGAGAACGTCCGG
>tag_2
AGGGCTAGGCGTGTCATTTGGCGAGGTCCTGAGGTCATGGAGCCA
>tag_3
CACTGGTCATCGACACTGGGAACCTGAGGTGAGTTCGCGCGCAAG
```  

`--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.


## Metrics

- Reads Mapped: R2 reads that successfully mapped to linker and tag-barcode.

- Reads Unmapped too Short: Unmapped R2 reads because read length < linker length + tag-barcode length.

- Reads Unmapped Invalid Linker: Unmapped R2 reads because of too many mismatches in linker sequence.

- Reads Unmapped Invalid Barcode: Unmapped R2 reads because of too many mismatches in tag-barcode sequence.
