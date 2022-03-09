## Features
- Map R2 reads to the tag barcode.

## Output
- raw_read_count.json: TODO
## Arguments
`--fq_pattern` R2 read pattern. The number after the letter represents the number of bases. The `fq_pattern` of CLindex is `L25C15`
`L` linker(common sequences)  
`C` tag barcode.

`--tag_barcode_fasta` Required. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < len(tag barcode) / 10 + 1. 
If no such tag exists, the read is classified as invalid.

`--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq` R2 read fastq.

