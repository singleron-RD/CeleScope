## Introduction

`multi_ffpe` is CeleScope's pipeline for single-cell FFPE RNA-Seq analysis. The workflow is highly similar to [single-cell RNA analysis](./multi_rna.md); only the assay-specific differences are described below. For the common usage, mapfile format, arguments, main output, and Seurat loading, please refer to [multi_rna.md](./multi_rna.md).

---

## Genome preparation

The main difference from `multi_rna` is in the genome preparation step — specifically, which `gene_biotype` entries should be kept in the GTF file.

For example, in mouse **Ensembl release 110**, the GTF file contains the following `gene_biotype` values:

```bash
awk -F'\t' '$3=="gene" {match($9, /gene_biotype "([^"]+)"/, arr); if(arr[1]!="") print arr[1]}' Mus_musculus.GRCm39.110.gtf | sort | uniq
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

When running `celescope utils mkgtf`, the following `gene_biotype` entries are kept by default:

```
gene_biotype=protein_coding,lncRNA,antisense,IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene;
```

Depending on the research needs, other `gene_biotype` categories can also be retained using the `--attributes` parameter. For example:

```bash
celescope utils mkgtf \
  --attributes "gene_biotype=protein_coding,lncRNA,antisense,IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene,miRNA,snoRNA,snRNA,scaRNA;" \
  Mus_musculus.GRCm39.110.gtf Mus_musculus.GRCm39.110.filtered.gtf
```

The remaining `mkref` steps are the same as for `multi_rna`.

---

## Usage

```bash
multi_ffpe \
    --mapfile ./ffpe.mapfile \
    --genomeDir {path to genomeDir} \
    --thread 16 \
    --soloCellFilter EmptyDrops_CR \
    --mod shell
```

The 4th column of the mapfile is optional — it is passed to `analysis` as the Celltypist model name for automatic cell-type annotation, same as `multi_rna`.

---

## miRNA analysis (optional)

The `mirna` step aligns the same FASTQ files against a **separate miRNA genome directory** and concatenates the miRNA counts with the gene expression matrix produced by `starsolo`, so that mRNA and miRNA features end up in the same `outs/raw` and `outs/filtered` matrices.

- It only runs when `--mirna_genomeDir` is provided; otherwise the step prints `Skip miRNA analysis.` and exits. Without it, `multi_ffpe` behaves exactly like `multi_rna`.
- STAR parameters are tuned for short reads: `--outFilterMismatchNoverLmax 0.1`, `--outFilterMatchNmin 15`, `--outFilterScoreMinOverLread 0.9`, `--outFilterMatchNminOverLread 0.9`, `--alignIntronMax 1`.

```bash
multi_ffpe \
    --mapfile ./ffpe.mapfile \
    --genomeDir {path to genomeDir} \
    --mirna_genomeDir {path to mirna genomeDir} \
    --thread 16 \
    --mod shell
```

---

## Other Built-in Modifications Different from the Standard scRNA Workflow

These modifications have already been applied in the code, so users don't need to change them manually.

- Changed the default value of `--outFilterMatchNmin` from 50 to 30.

- Added `SJ` to the default of `--soloFeatures` (`GeneFull_Ex50pAS Gene SJ`) to generate a splice junction matrix.

- Set the alignment strand to Reverse.