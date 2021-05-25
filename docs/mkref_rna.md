# mkref

## Features
- Create a genome reference directory.

## Input

- Genome fasta file

- Genome gtf file

- Mitochondria gene list file(optional)

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

## Parameters

`--genome_name` Required, genome name. 

`--fasta` Required, genome fasta file. 

`--gtf` Required, genome gtf file

`--mt_gene_list` Mitochondria gene list file. It is a plain text file with one gene per line. If not provided, will use `MT-` and `mt-` to determine mitochondria genes.

`--genomeDir` Default="./", output directory.

`--thread` Default=6, thread to use.

## Usage examples

- Human genome
```
celescope rna mkref \
 --genome_name Homo_sapiens_ensembl_99 \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.99.gtf
```

- Mouse genome
```
celescope rna mkref \
 --genome_name Mus_musculus_ensembl_99 \
 --fasta Mus_musculus.GRCm38.dna.primary_assembly.fa \
 --gtf Mus_musculus.GRCm38.99.gtf
```

