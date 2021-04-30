# cutadapt

## Features
- Trim adapters in R2 reads with cutadapt. Default adapters includes:
	- polyT=A{18}, 18 A bases. 
	- p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Illumina p5 adapter.

## Input
- R2 reads from step Barcode.

## Output
- `cutadapt.log` cutadapt output log file.

- `{sample}_clean_2.fq.gz` R2 reads file without adapters.

## Parameters

`--fq` R2 reads from step Barcode. Required. 

`--adapter_fasta` Addtional adapter Fasta file.

`--minimum_length` Default `20`. Discard processed reads that are shorter than LENGTH.

`--nextseq_trim` Default `20`. Quality trimming of reads using two-color chemistry (NextSeq). Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.

Since the regular quality-trimming algorithm cannot deal with this situation, you need to use the `--nextseq-trim` option: `cutadapt --nextseq-trim=20 -o out.fastq input.fastq`
This works like regular quality trimming (where one would use `-q 20` instead), except that the qualities of G bases are ignored.

`--overlap` Default `10`. Since Cutadapt allows partial matches between the read and the adapter sequence, short matches can occur by chance, leading to erroneously trimmed bases. For example, roughly 25% of all reads end with a base that is identical to the first base of the adapter. To reduce the number of falsely trimmed bases, the alignment algorithm requires that, by default, at least three bases match between adapter and read. This minimum overlap length can be changed globally (for all adapters) with the parameter `--overlap` (or its short version `-O`).

`--insert` Default `150`. Read2 insert length. To shorten each read down to a certain length, use the `--length` option or the short version `-l`.

`--gzip` Output gzipped fastq.

## Metrics
- Reads with Adapters : reads with sequencing adapters or reads two with poly A(read-through adpaters).

- Reads too Short : reads with read length less than 20bp after trimming.

- Reads Written : reads pass filtering.

- Base Pairs Processed : total processed base pairs.

- Base Pairs Quality-Trimmed : bases pairs removed from the end of the read whose quality is smaller than the given threshold.

- Base Pairs Written : base pairs pass filtering.