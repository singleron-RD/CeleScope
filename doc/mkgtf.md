The `celescope utils mkgtf` command is primarily used to refine genome annotation files by retaining only selected gene attributes.

For example, in the mouse Ensembl release 110 GTF file, the following `gene_biotype` values are present:

```bash
awk -F'\t' '$3=="gene" {match($9, /gene_biotype "([^"]+)"/, arr); if(arr[1]!="") print arr[1]}' Mus_musculus.GRCm39.110.gtf | sort | uniq
```

```text
IG_C_gene
IG_C_pseudogene
IG_D_gene
IG_D_pseudogene
IG_J_gene
IG_LV_gene
IG_pseudogene
IG_V_gene
IG_V_pseudogene
lncRNA
miRNA
misc_RNA
Mt_rRNA
Mt_tRNA
processed_pseudogene
protein_coding
pseudogene
ribozyme
rRNA
scaRNA
scRNA
snoRNA
snRNA
sRNA
TEC
transcribed_processed_pseudogene
transcribed_unitary_pseudogene
transcribed_unprocessed_pseudogene
translated_unprocessed_pseudogene
TR_C_gene
TR_D_gene
TR_J_gene
TR_J_pseudogene
TR_V_gene
TR_V_pseudogene
unitary_pseudogene
unprocessed_pseudogene
```

By default, the retained gene_biotype categories include protein-coding, lncRNA and multiple V(D)J-related gene segments.

Depending on the research objective, additional `gene_biotype` categories can be retained using the `--attributes` parameter. For example:

```bash
celescope utils mkgtf \
  --attributes "gene_biotype=protein_coding,lncRNA,IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene,miRNA,snoRNA,snRNA,scaRNA;" \
  Mus_musculus.GRCm39.110.gtf \
  Mus_musculus.GRCm39.110.filtered.gtf
```

If all `gene_biotype` categories should be retained, meaning no gene-level filtering is required, the `mkgtf` step can be skipped.

During the `celescope utils mkgtf` step, a default `mt_gene_list.txt` file is automatically generated. This file is created by identifying all gene symbols that begin with the prefix `"MT-"` or `"mt-"`, and is used as the default mitochondrial gene list in downstream analysis.

## Usage and Arguments

The command-line syntax is:

```bash
celescope utils mkgtf [options] <gtf> <out_gtf>
```

The main arguments are:

**`gtf`**
Path to the raw input GTF file. This argument is required.

**`out_gtf`**
Path to the processed output GTF file. This argument is required.

**`--attributes`**
A string specifying the gene attributes to retain.The default value is `--attributes "gene_biotype=protein_coding,lncRNA,antisense,IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene;"`

**`--skip_intron`**
By default, the tool adds intron entries to the output GTF based on the existing gene and exon structures. If this flag is provided, intron entries will not be added to the output GTF.
