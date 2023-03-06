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
```

```
$ cat celescope_genome.config
[genome]
genome_type = virus
fasta = EBV_genome.fasta
genome_name = EBV
```
## Arguments
`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exit.

`--fasta` Required. Genome fasta file. Use absolute path or relative path to `genomeDir`.

`--STAR_param` Additional parameters for the called software. Need to be enclosed in quotation marks. For example, `--{software}_param "--param1 value1 --param2 value2"`.

