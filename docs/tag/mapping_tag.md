## Features
- Align R2 reads to the tag barcode fasta.

## Output

- `{sample}_read_count.tsv` tab-delimited text file with 4 columns.

    `barcode` cell barcode  
    `tag_name`  tag name in barcode_fasta  
    `UMI`   UMI sequence  
    `read_count` read count per UMI  


## Arguments
`--fq_pattern` Required. R2 read pattern. The number after the letter represents the number of bases.         
`L` linker(common sequences)  
`C` tag barcode

`--barcode_fasta` Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < len(tag barcode) / 10 + 1. 
If no such tag exists, the read is classified as invalid.
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

`--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq` R2 read fastq.

