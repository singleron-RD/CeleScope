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

`--STAR_param` Additional parameters for the called software. Need to be enclosed in quotation marks. For example, `--{software}_param "--param1 value1 --param2 value2"`.

`--gtf` Required. Genome gtf file. Use absolute path or relative path to `genomeDir`.

`--mt_gene_list` Mitochondria gene list file. Use absolute path or relative path to `genomeDir`.
It is a plain text file with one gene per line. 
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.

`--genomeSAindexNbases` For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical 
value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal 
to 9, for 100 kiloBase genome, this is equal to 7.

