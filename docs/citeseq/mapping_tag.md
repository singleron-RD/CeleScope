## Features
- Assign tag to R2 reads.
## Output

- `{sample}_invalid_barcode.tsv` Reads count with invalid tag.

- `{sample}_read_count.tsv` Reads count with effective tag in each barcode.
## Arguments
`--fq_pattern` R2 read pattern. The number after the letter represents the number of bases. The `fq_pattern` of CLindex is `L25C15`
`L` linker(common sequences)  
`C` tag barcode.

`--barcode_fasta` Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < threshold. 
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

