## Features
- Create a fusion genome directory.

## Output

- STAR genome index files
- Genome config file

## Usage
```
celescope fusion mkref \
--genome_name {genome_name} \
--fasta fusion.fasta \
--fusion_pos fusion_pos.txt \
--genomeSAindexNbases 4
```


## Arguments
`--genomeDir` Default='./'. Output directory.

`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exit.

`--fasta` Fusion fasta file.

`--fusion_pos` fusion position file. A two column tab-delimited text file with header.
"pos" is the end postion of the first gene(1-based).
e.g.  
tag	pos  
PML_3	183  
PML_4	254  
PML_5	326  
PML_6	204

`--genomeSAindexNbases` STAR genomeSAindexNbases

