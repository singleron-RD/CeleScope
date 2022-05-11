## Features
- Create a virus genome reference directory.

## Output

- STAR genome index files
- Genome config file

## Usage
```
celescope capture_virus mkref \
    --genome_name EBV \
    --fasta EBV_genome.fasta \
    --genomeSAindexNbases 7
```

```
$ cat celescope_genome.config
[genome]
genome_type = virus
fasta = EBV_genome.fasta
genome_name = EBV
genomesaindexnbases = 7
```
## Arguments
`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exit.

`--fasta` Required. Genome fasta file. Use absolute path or relative path to `genomeDir`.

`--genomeSAindexNbases` STAR genomeSAindexNbases.

