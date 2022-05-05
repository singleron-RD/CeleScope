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

`--STAR_param` Additional parameters for the called software. Need to be enclosed in quotation marks. For example, `--{software}_param "--param1 value1 --param2 value2"`.

`--genomeSAindexNbases` For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical 
value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal 
to 9, for 100 kiloBase genome, this is equal to 7.

