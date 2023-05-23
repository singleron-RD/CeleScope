##Features
- Extract reads of mapping to VDJ gene.

For TCR library:
Total Reads: 10,000
Reads Mapped To Any V(D)J Genes: 9,000
Reads Mapped To TRA: 4,000
Reads Mapped To TRB: 5,000
For BCR library:
Total Reads: 20,000
Reads Mapped To Any V(D)J Genes: 15,000
Reads Mapped To IGH: 5,000
Reads Mapped To IGK: 5,000
Reads Mapped To IGL: 5,000

## Output
`bcrtcr_1.fq` barcode and umi of mapping to any V(D)J genes reads.
`bcrtcr_2.fq` sequence of mapping to any V(D)J genes reads.
## Arguments
`--fq2` raw fastq1 file path.

`--ref` reference name: hg38 or GRCm38.

`--seqtype` TCR or BCR.

`--out_dir` out directory.

