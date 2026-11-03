# CeleScope space (Spatial Transcriptomics)

## Overview

The `space` pipeline processes spatially barcoded transcriptome data generated with Singleron spatial kits. It handles read demultiplexing, barcode/UMI extraction, genome alignment via STARsolo, in-tissue spot filtering, and downstream analysis with Scanpy. The pipeline is invoked through `multi_space`, which generates and optionally runs one shell script per sample.

The pipeline supports three chemistries:

| Chemistry | Pattern (R1) | Description |
|---|---|---|
| `space-ff` | `U12-C8-L8-C8-L8-L16` | Fresh-frozen tissue; 96*96 spots |
| `space-ffpe` | `U12-C8-L8-C8-L8` | FFPE tissue, 96*96 spots |
| `space-ffpe-V1.1` | `U12-C8-L8-C8-L8` | FFPE tissue, version 1.1 144*144 spots |

Chemistry is auto-detected by default from the R1 read structure. For `space-ff`, the presence of the 16-nt linker 3 sequence `CTGTCTCTTATACACA` at positions 44–59 of R1 is the key signature.

---

## Prepare the genome directory

The spatial genome directory differs from the standard single-cell RNA-seq genome directory: the GTF annotation is **not filtered** by `celescope utils mkgtf`, so that short non-coding RNAs are retained, which is important for spatial QC.

### Homo sapiens

```bash
mkdir hs_ensembl_110
cd hs_ensembl_110

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz

conda activate celescope
celescope rna mkref \
 --genome_name Homo_sapiens_ensembl_110 \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.110.gtf
```

### Mus musculus

```bash
mkdir mmu_ensembl_110
cd mmu_ensembl_110

wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.110.gtf.gz

conda activate celescope

celescope rna mkref \
 --genome_name Mus_musculus_ensembl_110 \
 --fasta Mus_musculus.GRCm39.dna.primary_assembly.fa \
 --gtf Mus_musculus.GRCm39.110.gtf
```

---

## Prepare the spatial directory

The spatial directory contains the spot coordinate map, H&E image, and scale factors. It is generated from the imaging output using **[AtlasXbrowser](https://github.com/singleron-RD/AtlasXbrowser)**.

Required files inside the spatial directory:

| File | Description |
|---|---|
| `tissue_positions_list.csv` | Spot barcode, in-tissue flag (0/1), array row/column, pixel coordinates (no header) |
| `scalefactors_json.json` | Spot diameter, scale factors for image alignment |
| `figure/*.png` (or `.jpg`, `.jpeg`) | H&E tissue image |

The spatial directory is specified as the **4th column** of the mapfile. If the H&E image is missing, the pipeline will exit with an error.

---

## Prepare the mapfile

The mapfile is a **tab-delimited** text file with **four columns**. Each line represents a paired-end FASTQ pair (or a single pair for FFPE data, see Note below). Multiple lines with the same sample name (column 3) will be automatically concatenated before processing, which is useful when a sample was sequenced across multiple flow cells or lanes.

| Column | Description |
|---|---|
| 1 | Fastq file prefix (without `_1.fq.gz` / `_2.fq.gz` suffix) |
| 2 | Directory path containing the FASTQ files |
| 3 | Sample name (used as prefix for all output files) |
| 4 | Path to the spatial directory |

### Example: multiple FASTQ pairs per sample

Sample1 has 2 paired-end FASTQ pairs located in 2 different directories; Sample2 has 1 pair.

```bash
$cat ./sample.mapfile
fastq_prefix1	fastq_dir1	sample1	sample1_spatial_dir
fastq_prefix2	fastq_dir2	sample1	sample1_spatial_dir
fastq_prefix3	fastq_dir1	sample2	sample2_spatial_dir

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
```

> **Note for FFPE samples:** In FFPE data, the first 44 bp of R1 consist of barcode + linker + UMI, and the remaining sequence is cDNA used for analysis. R2 begins with a poly(T) stretch and the rest of R2 is typically low quality and not used. Therefore only the R1 FASTQ is required for analysis; R2 can be a placeholder or empty file.

---

## 4. Run the pipeline

Create a shell script `run.sh` in your working directory:

```bash
multi_space \
    --mapfile ./sample.mapfile \
    --genomeDir /path/to/hs_ensembl_110 \
    --thread 16 \
    --mod shell
```

After running `sh run.sh`, a `shell/` directory containing `{sample}.sh` files will be generated. Start the analysis:

```bash
bash ./shell/{sample}.sh
```

> **Important:** The `{sample}.sh` script must be executed from the **working directory** (the directory where `run.sh` was run), not from the `shell/` subdirectory.

### multi_space command-line options (key parameters)

| Option | Required | Default | Description |
|---|---|---|---|
| `--mapfile` | Yes | — | Path to the tab-delimited mapfile |
| `--genomeDir` | Yes | — | Path to the genome directory created by `celescope rna mkref` |
| `--thread` | No | 16 | Number of threads for STAR alignment. Using >20 threads is not recommended as STAR's mapping speed saturates above ~20 threads |
| `--mod` | No | `sjm` | Script type to generate: `sjm` (Simple Job Manager) or `shell` |


---

## Pipeline steps

The `space` pipeline consists of four sequential steps:

### Step 1 — `sample`

Performs sample concatenation (merges multiple FASTQ pairs for the same sample) and basic file existence checks.

### Step 2 — `starsolo`

The core alignment and quantification step. Key operations:

1. **Read trimming**: A 12-nt poly(A) 3′ adapter (`AAAAAAAAAAAA`) is clipped from the 3′ end of reads.
2. **STAR alignment**: Reads are aligned to the reference genome with STAR 2.7.11a in **STARsolo** mode.
3. **Barcode/UMI extraction**: Two spot-barcode segments (positions 0–7 and 16–23 in R1, for `space-ff`) are matched against pre-defined whitelists allowing up to two mismatches in total per barcode (`soloCBmatchWLtype EditDist_2`). The UMI spans positions 32–43.
4. **Filtering disabled**: Automatic cell/spot calling is disabled (`--soloCellFilter None`) because spatial arrays have a fixed, known number of spots and every spot is expected to have some signal.
5. **Multi-feature counting**: Reads are counted simultaneously at two feature levels:
   - `GeneFull_Ex50pAS`: pre-mRNA-permissive counts (a read is assigned to a gene if ≥50% of it overlaps an exon; intronic reads are also counted)
   - `Gene`: canonical exon-only counts
   The `GeneFull_Ex50pAS` counts are used as the primary matrix for all downstream analysis and reporting.
6. **In-tissue filtering**: Only spots annotated as `in_tissue = 1` in the spatial manifest are retained from the raw matrix.
7. **Output**: raw matrix, in-tissue filtered matrix, sorted BAM, and a QC metrics HTML report.

> **For FFPE chemistry**, the trimming differs: `--soloStrand Reverse --clip5pNbases 44` is applied instead, and only R1 is used as the genomic read (R2 is single-end placeholder).

#### starsolo options (key parameters)


Other less commonly used options: `--soloCBmatchWLtype` (barcode mismatch tolerance, default `EditDist_2`), `--report_soloFeature` (feature type used in the report, default `GeneFull_Ex50pAS`), `--outSAMtype` (default `BAM SortedByCoordinate`, set `None` to skip BAM), `--STAR_param` (additional raw STAR arguments), `--adapter_3p` (3′ adapter, default `AAAAAAAAAAAA`). `--pattern` and `--whitelist` are only needed with `--chemistry customized`.

### Step 3 — `mirna` (optional)

Performs a separate miRNA-level alignment and quantification using a dedicated miRNA genome index. This step is skipped by default. To enable it, prepare a miRNA genome directory (e.g. with `celescope rna mkref` using a miRNA fasta) and add `--mirna_genomeDir /path/to/mirna_genome` to the `mirna` step in the shell script. miRNA analysis uses `--outFilterMatchNmin 15` and `--soloCellFilter None`.

### Step 4 — `analysis`

Downstream spatial analysis using Scanpy. Key operations:

1. **Matrix conversion**: The raw and in-tissue-filtered MTX matrices are converted to 10x Visium-compatible HDF5 files (`raw_feature_bc_matrix.h5` and `filtered_feature_bc_matrix.h5`).
2. **Spatial directory generation**: The input spatial manifest is copied to `outs/spatial/`, with the H&E image resized to 2,000 px (hires) and 600 px (lowres), scale factors augmented with `microns_per_pixel` and `regist_target_img_scalef`, and the coordinate CSV renamed for Visium compatibility.
3. **Data loading**: The raw matrix is loaded with `sc.read_visium`.
4. **Filtering**: Genes detected in fewer than 3 spots, and spots with fewer than 1 gene (configurable via `--min_genes`), are removed.
5. **QC metrics**: Mitochondrial gene fraction is computed and reported per-spot.
6. **Normalization (Pearson residuals)**: Analytic Pearson residuals are computed (`scanpy.experimental.pp.normalize_pearson_residuals`), providing a variance-stabilizing normalization analogous to Seurat's SCTransform.
7. **Highly variable genes**: 2,000 highly variable genes are selected with the `pearson_residuals` flavor.
8. **Dimensionality reduction**: PCA (50 components), k-nearest-neighbor graph (15 neighbors, 25 PCs), and UMAP embedding are computed.
9. **Clustering**: Leiden community detection at resolution 0.3.
10. **Differential expression**: Cluster marker genes are identified by Wilcoxon rank-sum test (adjusted *P* < 0.05) using size-factor-normalized (10,000 counts/spot) log1p-transformed counts.


---

## Loading results into Seurat and Scanpy

The `filtered_feature_bc_matrix.h5` and `spatial/` folder follow the [10x Space Ranger format](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs) and can be loaded directly into Seurat.

The examples below use `./outs` as the data directory, following the structure created by the copy commands below. The Visium format requires `filtered_feature_bc_matrix.h5` and the `spatial/` folder to be siblings, so gather them into one directory first:

```bash
mkdir -p ./outs
cp -r ./outs/{sample}/analysis/spatial ./outs/spatial
cp ./outs/{sample}/analysis/filtered_feature_bc_matrix.h5 ./outs/
```

### Seurat v5

```R
library(Seurat) # v5
obj <- Load10X_Spatial(data.dir = "./outs")
obj@meta.data$orig.ident = "sample1"
```

### Seurat v4

Seurat v4 does not recognise the `tissue_positions.parquet` file in the spatial directory and requires the `tissue_positions_list.csv` file. Rename the files before loading:

```bash
# rename parquet to avoid conflict
mv outs/spatial/tissue_positions.parquet outs/spatial/parquet.bak
mv outs/spatial/positions_list.csv outs/spatial/tissue_positions_list.csv
```

```R
library(Seurat) # v4
obj <- Load10X_Spatial(data.dir = "./outs")
obj@meta.data$orig.ident = "sample1"
```

### Scanpy

```python
import scanpy as sc

adata = sc.read_h5ad("./outs/{sample}/analysis/rna.h5ad")
```

The `.h5ad` file produced by the `analysis` step (`outs/{sample}/analysis/rna.h5ad`) contains:

- `adata.X`: Pearson-residual-normalized expression (used for dimensionality reduction)
- `adata.layers["counts"]`: raw counts
- `adata.layers["normalized"]`: size-factor-normalized (10,000 counts/spot) log1p-transformed expression (used for DE testing)
- `adata.obs["cluster"]`: Leiden cluster assignments
- `adata.obs["total_counts"]`: total UMI per spot
- `adata.obs["pct_counts_mito"]`: mitochondrial UMI fraction
- `adata.obsm["X_pca"]`, `adata.obsm["X_umap"]`: PCA and UMAP embeddings
- Spatial coordinates and image data are preserved in `adata.uns["spatial"]`


