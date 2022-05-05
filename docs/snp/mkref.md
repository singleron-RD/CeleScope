## Features
- Create dictionary file and fasta index for gatk SplitNCigarReads.
(https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) 
Need to run `celescope rna mkref` first

## Output
- fasta index
- gatk dictionary file

## Usage
```
# run celescope rna mkref first
celescope snp mkref \
 --genome_name Homo_sapiens_ensembl_99 \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
```
## Arguments
`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exit.

`--fasta` Required. Genome fasta file. Use absolute path or relative path to `genomeDir`.

`--STAR_param` Additional parameters for the called software. Need to be enclosed in quotation marks. For example, `--{software}_param "--param1 value1 --param2 value2"`.

