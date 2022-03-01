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
`C` tag barcode.

`--barcode_fasta` Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < len(tag barcode) / 10 + 1. 
If no such tag exists, the read is classified as invalid.

You can find the barcode fasta file under `celescope/data/Clindex`
```
>CLindex_TAG_1
CGTGTTAGGGCCGAT
>CLindex_TAG_2
GAGTGGTTGCGCCAT
>CLindex_TAG_3
AAGTTGCCAAGGGCC
>CLindex_TAG_4
TAAGAGCCCGGCAAG
>CLindex_TAG_5
TGACCTGCTTCACGC
>CLindex_TAG_6
GAGACCCGTGGAATC
>CLindex_TAG_7
GTTATGCGACCGCGA
>CLindex_TAG_8
ATACGCAGGGTCCGA
>CLindex_TAG_9
AGCGGCATTTGGGAC
>CLindex_TAG_10
TCGCCAGCCAAGTCT
>CLindex_TAG_11
ACCAATGGCGCATGG
>CLindex_TAG_12
TCCTCCTAGCAACCC
>CLindex_TAG_13
GGCCGATACTTCAGC
>CLindex_TAG_14
CCGTTCGACTTGGTG
>CLindex_TAG_15
CGCAAGACACTCCAC
>CLindex_TAG_16
CTGCAACAAGGTCGC
```

`--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq` R2 read fastq.

