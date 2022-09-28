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

You can find the barcode fasta file under `celescope/data/sweetseq`.

`--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fq` R2 read fastq.

