## Features
- Align R2 reads to the tag barcode fasta.

## Output

- `{sample}_read_count.tsv` tab-delimited text file with 4 columns.

    `barcode` cell barcode  
    `tag_name`  tag name in barcode_fasta  
    `UMI`   UMI sequence  
    `read_count` read count per UMI

- `{sample}_invalid_barcode.tsv` tab-delimited text file with 2 columns.
    `tag_barcode` tag barcodes that do not match with any sequence in `--barcode_fasta`.
    `read_count` invalid tag barcode read counts
## Arguments
`--fq_pattern` R2 read pattern. The number after the letter represents the number of bases. CLindex is `L25C15` and sweetseq is `L23C15`
`L` linker(common sequences)  
`C` tag barcode.

`--barcode_fasta` Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < threshold. 
If no such tag exists, the read is classified as invalid.

You can find the barcode fasta file under `celescope/data/sweetseq`.

`--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq` R2 read fastq.

