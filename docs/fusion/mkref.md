

## Arguments
`--genomeDir` Default='./'. Output directory.

`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exist.

`--fasta` fusion fasta file

`--fusion_pos` fusion position file. A two column tab-delimited text file with header.
"pos" is the end postion of the first gene(1-based).
e.g.  
tag	pos  
PML_3	183  
PML_4	254  
PML_5	326  
PML_6	204

`--genomeSAindexNbases` STAR genomeSAindexNbases

