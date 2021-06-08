## Features
- Create a genome reference directory.

## ## Output

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
`--genomeDir` Default='./'. Output directory.

`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exist.

`--fasta` Required. Genome fasta file.

`--gtf` Required. Genome gtf file.

`--mt_gene_list` Mitochondria gene list file. It is a plain text file with one gene per line. 
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.

