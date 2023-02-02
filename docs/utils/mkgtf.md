##Features

- Filter GTF files.

GTF files downloaded from sites like ENSEMBL and UCSC often contain transcripts and genes which need to be filtered.
Celescope provides mkgtf, a simple utility to filter genes. The command syntax requires input and output fitlered GTF file names.

The following gene biotypes will be kept.
```
protein_coding 
lncRNA 
antisense 
IG_LV_gene 
IG_V_gene 
IG_V_pseudogene 
IG_D_gene 
IG_J_gene 
IG_J_pseudogene 
IG_C_gene 
IG_C_pseudogene 
TR_V_gene 
TR_V_pseudogene 
TR_D_gene 
TR_J_gene 
TR_J_pseudogene 
TR_C_gene
```
## Arguments
`--attributes` Attributes to keep. Example: `gene_biotype=protein_coding,lncRNA,antisense;`.

