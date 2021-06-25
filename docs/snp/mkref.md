## Features
- Create dictionary file and fasta index for gatk SplitNCigarReads.
(https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) 
Need to run `celescope rna mkref` first

## Output
- fasta index
- gatk dictionary file


## Arguments
`--genomeDir` Default='./'. Output directory.

`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exit.

`--fasta` fasta file

