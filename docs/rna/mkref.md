## Features
- Create a genome reference directory.

## Output

- STAR genome index files

- Genome refFlat file

- Genome config file
```
$ cat celescope_genome.config
[genome]
genome_name = Homo_sapiens_ensembl_99
genome_type = rna
fasta = Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf = Homo_sapiens.GRCh38.99.gtf
refflat = Homo_sapiens_ensembl_99.refFlat
```
## Arguments
`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exit.

`--fasta` Required. Genome fasta file. Use absolute path or relative path to `genomeDir`.

`--gtf` Required. Genome gtf file. Use absolute path or relative path to `genomeDir`.

`--mt_gene_list` Mitochondria gene list file. Use absolute path or relative path to `genomeDir`.
It is a plain text file with one gene per line. 
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.

`--genomeSAindexNbases` STAR genomeSAindexNbases.

