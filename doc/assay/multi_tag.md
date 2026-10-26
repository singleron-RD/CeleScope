## Introduction to CLindex®

[CLindex® Sample Multiplexing Kit](https://singleron.bio/products/clindex/) is a sample multiplexing kit developed by Singleron Biotechnologies for single-cell RNA-seq experiments. Prior to single-cell capture, a known DNA tag sequence is covalently attached to all cells in each sample via a chemical reaction. The tagged cells are then pooled into a single GEXSCOPE® library construction reaction, where the tag sequences are co-sequenced together with the cellular transcriptome.

During data analysis, the CeleScope `multi_tag` pipeline extracts tag sequences from R2 reads and aligns them against a reference FASTA to determine the sample (or treatment group) of origin for each cell. In this way, up to **16 samples** can be processed in a single run with the 16-plex Kit, and higher-order encoding can be achieved by combining multiple orthogonal tags. This dramatically increases cell throughput while reducing batch effects and sequencing cost per sample.

The `multi_tag` pipeline provides two tag-assignment algorithms (**SNR** and **Otsu**), integrates tag assignments with scRNA-Seq clustering results, and can split downstream files — expression matrices, FASTQ, BAM, and VDJ tables — by tag.

## Beyond CLindex® — general-purpose cell identity assignment

Although originally designed for CLindex® sample multiplexing, `multi_tag` is a **general-purpose pipeline for assigning cell identities based on any kind of DNA tag sequence that is co-sequenced alongside the cellular transcriptome**. The pipeline operates purely on the relationship between R2 read sequences and a user-supplied tag barcode FASTA reference, and does not make any assumptions about the biological origin of those tags. As a result, it can be directly applied — with only a custom `--barcode_fasta` and `--fq_pattern` — to a wide range of experiments that require per-cell identity calling, including but not limited to Perturb-seq / CRISPR.

For these use cases, simply replace the CLindex®-specific `--barcode_fasta` with your own tag reference FASTA, and adjust `--fq_pattern` to match the structure of your R2 reads. The rest of the pipeline — tag alignment, UMI collapsing, SNR/Otsu cell calling, visualization, and data splitting — works identically.

---

## Usage

### Prerequisites

Before running `multi_tag`, you must first complete the scRNA-Seq analysis using `multi_rna`. The tag pipeline depends on the cell barcode calling and clustering results from the RNA pipeline.

---

### Generate scripts for each sample

In your working directory, create a shell script named `run.sh`:

```bash
multi_tag \
  --mapfile ./tag.mapfile \
  --barcode_fasta ./CLindex_TAG.fasta \
  --fq_pattern L25C15 \
  --mod shell
```

To also generate per-tag split expression matrices, add `--split_matrix`:

```bash
multi_tag \
  --mapfile ./tag.mapfile \
  --barcode_fasta ./CLindex_TAG.fasta \
  --fq_pattern L25C15 \
  --split_matrix \
  --mod shell
```

#### Arguments

`--mapfile`
Required. The mapfile is a tab-delimited text file containing at least 4 columns. Each line represents one pair of FASTQ files.

* 1st column: FASTQ file prefix
* 2nd column: FASTQ file directory
* 3rd column: Sample name, which will be used as the prefix for all output files
* 4th column: Required, matched directory (`match_dir`). The output directory path from the previously completed CeleScope scRNA-Seq run.

Example:

Sample 1 has two paired-end FASTQ files located in two different directories, `fastq_dir1` and `fastq_dir2`. Sample 2 has one paired-end FASTQ file located in `fastq_dir1`.

```bash
$ cat ./tag.mapfile
fastq_prefix1	fastq_dir1	sample1	/path/to/rna/sample1
fastq_prefix2	fastq_dir2	sample1	/path/to/rna/sample1
fastq_prefix3	fastq_dir1	sample2	/path/to/rna/sample2

$ ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
```

`--barcode_fasta`
Required. Path to a FASTA file containing the expected tag barcode sequences. The pipeline checks mismatches between tag barcode sequences in R2 reads and all sequences in this file. A read is assigned to the tag with mismatch < threshold. If no such tag exists, the read is classified as invalid.

You can find example barcode FASTA files under `celescope/data/Clindex` or `celescope/data/sweetseq`.

Example content:

```
>CLindex_TAG_1
CGTGTTAGGGCCGAT
>CLindex_TAG_2
GAGTGGTTGCGCCAT
>CLindex_TAG_3
AAGTTGCCAAGGGCC
```

`--fq_pattern`
The structural pattern of the Read 2 sequence. The number after each letter represents the number of bases.

* `L` — Linker (common sequence)
* `C` — Tag barcode

For CLindex™ kits, the default is `L25C15` (Linker 25 bp, Barcode 15 bp).

`--mod`
Specifies the script type to generate. Available options include `sjm`, which uses [Simple Job Manager](https://github.com/StanfordBioinformatics/SJM), and `shell`, which generates standard shell scripts.

After running:

```bash
sh run.sh
```

a `shell` directory containing `{sample}.sh` files will be generated.

Start the analysis by running:

```bash
bash ./shell/{sample}.sh
```

Please note that `./shell/{sample}.sh` must be executed from the working directory. It should not be run from inside the `shell` directory.

---

## Pipeline Steps

### Step 1. barcode

Demultiplex R1 reads: extract cell barcodes and UMIs from R1, correct barcodes against a whitelist, and write corrected R2 with demultiplexed read names.

#### Process

Detect single-cell chemistry (default `auto`) and parse the R1 pattern into cell barcode (`C`), linker (`L`), UMI (`U`), and poly-T (`T`) segments against the corresponding whitelist. Build a precomputed lookup table allowing up to 1 Hamming distance from each whitelist barcode so that observed barcode sequences can be corrected in constant time. For each paired R1–R2 read, extract the cell barcode and UMI from R1, correct the barcode via the lookup table, and write the read to a demultiplexed R2 FASTQ whose header is reformatted as `corrected_barcode:umi:read_index`. Reads whose barcode cannot be corrected are discarded.

#### Arguments

`--chemistry`
Default: `auto`. Single-cell chemistry. Choices: keys from `chemistry_dict` (e.g. `auto`, `GEXSCOPE-V1`, `GEXSCOPE-V2`, `GEXSCOPE-V3`, `flv`, etc.).

`--pattern`
Optional. R1 read pattern (e.g. `C8L16C8L16C8L1U12T18`). The number after each letter represents the number of bases. Only needed when `--chemistry` is `auto` or overriding the default.
* `C` — cell barcode
* `L` — linker (common sequence)
* `U` — UMI
* `T` — poly-T

`--whitelist`
Optional. Path to a cell barcode whitelist file, one barcode per line. Only needed when overriding the chemistry's built-in whitelist.

#### Output

`{outdir}/{sample}_2.fq`
Demultiplexed R2 FASTQ. Only reads with valid (or correctable) cell barcodes are kept. The read header is reformatted as `{corrected_barcode}:{umi}:{read_index}`. This file is consumed by `mapping_tag`.

#### Metrics (added to HTML report)

| Metric Name | Description |
|---|---|
| Raw Reads | Total paired reads input |
| Valid Reads | Reads whose cell barcode exactly matches the whitelist |
| Corrected Reads | Reads whose cell barcode was not in the whitelist but was within 1 Hamming distance and was corrected |
| Q30 of Barcode | Fraction of cell barcode bases with Phred quality score ≥ 30 |
| Q30 of UMI | Fraction of UMI bases with Phred quality score ≥ 30 |

---

### Step 2. mapping_tag

Map R2 reads to the tag barcode FASTA reference.

#### Process

Parse `fq_pattern` to extract linker (`L`) and tag barcode (`C`) segments from the R2 read sequence. Determine a per-segment mismatch threshold, either specified on the command line or auto-calculated from the size of each FASTA reference (≤ 2 mismatches, capped at the segment length). Build a precomputed lookup table per segment that maps every possible sequence within the threshold to its originating tag, enabling O(1) lookup for every R2 read. When `--linker_fasta` is provided, reads whose linker fails a Hamming-distance check (`floor(len(linker) / 10) + 1`) are discarded. For the remaining reads, extract the tag barcode segments and look each up in its corresponding table; a read is assigned to the tag formed by joining per-segment matches with underscores, and its UMI count is incremented accordingly. Reads with one or more non-matching segments are classified as invalid.

#### Arguments

`--fq_pattern`
Default: `L25C15`. R2 read pattern. The number after each letter represents the number of bases.
* `L` — linker (common sequence)
* `C` — tag barcode

`--barcode_fasta`
Required. Tag barcode FASTA file. Comma-separated if there are multiple `C` segments in `fq_pattern`. The number of FASTA files must equal the number of `C` segments.

`--linker_fasta`
Optional. If provided, it checks mismatches between the linker sequence in R2 reads and all linker sequences in `linker_fasta`. If no match is found within the Hamming distance threshold, the read is classified as invalid.

`--mismatch`
Optional. Maximum number of mismatches allowed between tag barcode in R2 reads and `barcode_fasta`. Comma-separated if there are multiple `C` segments (one value per segment). If a single value is given, it applies to all `C` segments. If not specified, auto-determined per segment based on barcode count:
```
> 100,000 barcodes: 0
10,001 - 100,000 barcodes: 1
≤ 10,000 barcodes: 2
```

#### Output

`{outdir}/{sample}_read_count.tsv`
Tab-delimited text file with 4 columns.

| Column | Description |
|---|---|
| `barcode` | Cell barcode extracted from the demultiplexed R2 read name |
| `tag_name` | Assigned tag name. When multiple `C` segments exist, per-segment names are joined with `_` |
| `UMI` | UMI sequence extracted from the demultiplexed R2 read name |
| `read_count` | Number of reads for this (barcode, tag_name, UMI) combination |

`{outdir}/{sample}_invalid_barcode.tsv`
Tab-delimited text file with invalid barcode information. The header is dynamically generated based on the number of `C` segments. For each segment `i` there are two columns:

| Column | Description |
|---|---|
| `tag_barcode_1` ... `tag_barcode_N` | Observed tag barcode sequence for each `C` segment |
| `valid_1` ... `valid_N` | `"True"` if the sequence matched the reference (within mismatch threshold), `"False"` otherwise |
| `read_count` | Number of reads with this invalid barcode tuple |

Rows are sorted by `read_count` descending.

#### Metrics (added to HTML report)

| Metric Name | Description |
|---|---|
| Reads Mapped | R2 reads that successfully mapped to both linker (if applicable) and tag barcode |
| Reads Unmapped Invalid Linker | Reads discarded due to linker mismatch |
| Reads Unmapped Invalid Barcode | Reads discarded because tag barcode failed to match any reference |

---

### Step 3. count_tag

Assign tag identity to each cell barcode and summarize cell counts per tag.

#### Process

Filter mapped reads to scRNA-Seq filtered cell barcodes, then collapse reads by cell and tag to obtain per-cell per-tag UMI counts. Every scRNA-Seq cell is included even if it has no tag UMIs. Two tag-assignment methods are available.

The **SNR method** is the default and works well for most experiments. It relies on the relative gap between the top and the second tag within each cell, so it is robust when all tag intensities are comparable. When tags show highly uneven labeling strength — for example, one tag consistently produces many more UMIs than the others — the SNR of cells carrying the weaker tag may become indistinguishable from background. In that case, consider switching to the **Otsu method**, which fits an independent threshold per tag and therefore tolerates large absolute differences in tag intensity.

**SNR method (`--tag_threshold_method snr`, default)**

**1. Per-cell signal-to-noise ratio (SNR)**

For each cell, let its per-tag UMI counts be sorted in descending order: x₁ ≥ x₂ ≥ ... ≥ xₖ. With `dim = d` (default 1),

```
             ⎧ x_d / x_(d+1)     if x_d > 0 and x_(d+1) > 0
  SNR =      ⎨ ∞                 if x_(d+1) = 0  (no noise)
             ⎩ 0                 if x_d = 0
```

The signal is the d-th highest tag and the noise is the (d+1)-th highest. With `dim=1`, SNR is simply the ratio of the top tag to the second tag.

**2. Automatic `UMI_min`**

```
UMI_min = max( min( percentile₅(UMIs), median(UMIs) / 10 ), 1 )
```

where UMIs are the per-cell total tag UMI counts across all scRNA-Seq cells. A cell whose total tag UMI is below this threshold is too sparse for reliable calling.

**3. Automatic `snr_min`**

Restrict to cells satisfying UMI ≥ UMI_min. If the median SNR among these cells is infinite (no noise observed), set snr_min = 10; otherwise,

```
snr_min = max( median(SNRs) × coefficient, 2 )
```

The default `coefficient` is 0.1; a smaller value produces fewer Multiplet calls.

**4. Cell classification**

| Condition | Assignment |
|---|---|
| UMI < UMI_min | **Undetermined** — too few tag UMIs |
| UMI ≥ UMI_min and SNR < snr_min | **Multiplet** — multiple tags at comparable levels |
| UMI ≥ UMI_min and SNR ≥ snr_min | Assigned to the top `dim` tag names (joined with `_` when `dim > 1`) |

**Otsu method (`--tag_threshold_method otsu`)**

Compute a positive/negative threshold for each tag independently, then classify each cell by counting how many of its tags are above their respective thresholds.

**1. Per-tag threshold via Otsu's method**

For each tag's UMI-count column:
* Remove cells with zero UMI for that tag; if fewer than 50 cells remain, set the threshold to 1.
* Log10-transform the remaining counts.
* Build a histogram of the log-space values with bin width 0.2, spanning from 0 to max(log-counts) + 0.2.
* Find the bin threshold τ that maximizes the between-class variance:

```
w₁(τ) × w₂(τ) × (μ₁(τ) − μ₂(τ))²
```

where w₁, μ₁ are the weight and mean of counts falling below τ, and w₂, μ₂ are those above τ.
* The final threshold in UMI space is ⌈10^τ⌉, i.e. the ceiling of 10 raised to the log-space threshold.

**2. Cell classification**

For each cell, compare each tag's UMI count to its per-tag Otsu threshold:

| Tags above their threshold | Assignment |
|---|---|
| 0 | **Undetermined** — no tag is significantly expressed |
| ≥ 2 | **Multiplet** — multiple tags are significantly expressed |
| exactly 1 | Assigned to that tag |

Both methods conclude by summarizing the number and percentage of cells falling into `Undetermined`, `Multiplet`, and each specific tag.

#### Arguments


`--UMI_min`
Default: `auto`. Minimum UMI threshold. Cell barcodes with valid UMI below this are classified as *Undetermined*. When `auto`, it is derived from the distribution of per-cell total tag UMI.

`--snr_min`
Default: `auto`. Minimum signal-to-noise ratio. Cell barcodes with UMI ≥ `UMI_min` but SNR below this are classified as *Multiplet*. When `auto`, it is derived from the distribution of SNR values.

`--coefficient`
Default: `0.1`. Only used when `snr_min` is `auto`. A smaller coefficient results in fewer *Multiplet* assignments.

`--dim`
Default: `1`. Number of expected orthogonal tag dimensions per cell.

`dim` controls two behaviors in the **SNR method** (it has **no effect** on the Otsu method):

1. **SNR calculation**. With per-cell tag UMI counts sorted descending as x₁ ≥ x₂ ≥ ... ≥ xₖ, SNR is computed as x_d / x_(d+1) where d = `dim`. The signal is the d-th highest tag and the noise is the (d+1)-th. With `dim=1` (default), SNR is the ratio of highest tag to the second-highest tag — optimal when each cell should carry exactly **one** tag. With `dim=2`, SNR becomes the ratio of the **second** tag to the third — appropriate when each cell legitimately carries **two** orthogonal tags and you want the pipeline to distinguish true double-positive cells from background (e.g., combinatorial 2D CLindex, dual-pool hashing, or Perturb-seq with two orthogonal gRNA sets).

2. **Final tag assignment**. When a cell passes the SNR threshold, its identity is the top d tag names joined with `_`. For example with `dim=2`, a cell whose two highest tags are `A` and `B` is labeled `A_B` rather than just `A`.

Set `dim` to the number of independent tag encoding dimensions in your experiment. For CLindex® single-plex or simple Perturb-seq / hashing experiments, keep the default `dim=1`.

> **Note — Otsu method limitation**: `dim` is **only honored by the SNR method**. The Otsu method currently treats any cell with ≥2 tags above their respective per-tag thresholds as a *Multiplet* and cannot join multiple positive tag names into a combined identity (e.g. `A_B`). Therefore Otsu is only appropriate for experiments where `dim=1` (one tag per cell). If your experiment has `dim≥2`, you must use the SNR method.

`--tag_threshold_method`
Default: `snr`. Method to determine tag positivity for each cell. Choices:
* `snr` — use SNR and `UMI_min` (original method).
* `otsu` — use Otsu's method to calculate a threshold for each tag independently. Cells with multiple tags above their respective thresholds are classified as *Multiplet*. Cells with no tags above threshold are classified as *Undetermined*.


#### Output

`{outdir}/{sample}_umi_tag.tsv`
Tab-delimited file with one row per scRNA-Seq cell barcode (index). Columns:
* First column (index): `barcode`
* Middle columns: UMI count per tag (column names match tag names in `barcode_fasta`)
* Last column (`tag`): Assigned tag — one of `Undetermined`, `Multiplet`, or a specific tag name

`{outdir}/otsu_plot/` (only when `--tag_threshold_method otsu`)

| File | Description |
|---|---|
| `{tag_name}_otsu.png` | Histogram of log10-transformed UMI counts for the tag, with the Otsu threshold marked as a red vertical line |
| `threshold.tsv` | Two columns: `tag` (tag name) and `threshold` (the Otsu-calculated threshold as an integer) |

#### Metrics (added to HTML report)

| Metric Name | Description |
|---|---|
| Mapped Reads in Cells | Mapped reads whose cell barcode is in the scRNA-Seq filtered cell set |
| Median UMI per Cell | Median total tag UMI across all scRNA-Seq cell barcodes |
| Mean UMI per Cell | Mean total tag UMI across all scRNA-Seq cell barcodes |
| Tag Cells (table) | For `Undetermined`, `Multiplet`, and each tag: cell count and percentage of total scRNA-Seq cells |

---

### Step 4. analysis_tag

Combine scRNA-Seq clustering information with tag assignment and generate t-SNE visualizations.

#### Process

Read t-SNE coordinates and cluster labels from the matched scRNA-Seq run, merge with per-cell tag assignment on cell barcode, and generate t-SNE scatter plots colored by cluster and by tag. Write the merged table to `{sample}_tsne_tag.tsv`.


#### Output

`{outdir}/{sample}_tsne_tag.tsv`
Tab-delimited file with one row per cell barcode. It contains all columns from the scRNA-Seq t-SNE file (typically `tSNE_1`, `tSNE_2`, `cluster`, `Gene_Counts`) plus one UMI-count column per tag and the final `tag` assignment column.

#### HTML Report Plots

| Plot | Description |
|---|---|
| t-SNE plot colored by cluster | Plotly scatter plot of t-SNE coordinates colored by scRNA-Seq Leiden cluster |
| t-SNE plot colored by tag | Plotly scatter plot of t-SNE coordinates colored by tag assignment (`Undetermined`, `Multiplet`, or specific tag name) |

---

### Step 5. split_tag

Split scRNA-Seq data according to tag assignment.

#### Process

Group cell barcodes by their assigned tag. For every tag, optionally split each requested data type: scRNA-Seq expression matrices (`--split_matrix`), FASTQ files (`--split_fastq`), aligned BAM files (`--split_bam`), VDJ count tables (`--split_vdj`), and full-length VDJ annotations (`--split_fl_vdj`). All split files follow the format of their input counterparts (e.g. 10x Genomics matrix format, coordinate-sorted BAM, VDJ confident-cell and clonotype tables). `Undetermined` and `Multiplet` are treated as regular tags and will also be split if present in the data.

#### Arguments

`--split_matrix`
If set, split scRNA-Seq matrix file according to tag assignment.

`--split_fastq`
If set, split scRNA-Seq FASTQ files according to tag assignment.

`--split_bam`
If set, split scRNA-Seq featureCounts BAM file according to tag assignment. Use together with `--bam_file`.

`--bam_file`
Path to the scRNA-Seq BAM file to split. Required only when `--split_bam` is used.

`--split_vdj`
If set, split scRNA-Seq VDJ count file according to tag assignment.

`--split_fl_vdj`
If set, split scRNA-Seq full-length VDJ annotation, FASTA, and clonotypes file according to tag assignment.

`--vdj_dir`
Match CeleScope VDJ directory. Required when `--split_vdj` or `--split_fl_vdj` is specified.

`--match_dir`
Match CeleScope scRNA-Seq directory. Required for `--split_fastq`. Also used as the fallback for matrix location when `--matrix_dir` is not provided.

`--matrix_dir`
Match CeleScope scRNA-Seq matrix directory. Alternative to obtaining the matrix via `--match_dir`.

`--R1_read`
Comma-separated R1 read paths. Required when `--split_fastq` is specified.

#### Output

`{outdir}/matrix/{tag}_filtered_feature_bc_matrix/` (when `--split_matrix`)
Per-tag expression matrices in 10x Genomics format. Subfiles:
* `matrix.mtx.gz` — Sparse expression matrix (barcodes × genes)
* `barcodes.tsv.gz` — Cell barcodes for this tag
* `features.tsv.gz` — Gene IDs and gene symbols

`{outdir}/fastqs/` (when `--split_fastq`)
Per-tag FASTQ files:
* `{tag}_1.fq` — R1 reads
* `{tag}_2.fq` — R2 reads

`{outdir}/bam/{tag}.bam` (when `--split_bam`)
Coordinate-sorted BAM file containing reads assigned to each tag. The header is copied from the input BAM; unmapped reads are skipped.

`{outdir}/vdj/` (when `--split_vdj`)
Per-tag VDJ count files. The `{seqtype}` is determined from the chain column (`TCR` for `TRA/TRB`, `BCR` for `IGH/IGL/IGK`).
* `{tag}_{seqtype}_cell_confident.tsv` — Filtered VDJ confident-cell table for this tag
* `{tag}_{seqtype}_clonotypes.tsv` — Clonotype table with columns `clonetype_ID`, chain-specific `aaSeqCDR3` and `nSeqCDR3`, `barcode_count`, and `percent`

`{outdir}/fl_vdj/` (when `--split_fl_vdj`)
Per-tag full-length VDJ files.
* `{tag}_{seqtype}_contig_annotations.csv` — Filtered annotation CSV for this tag
* `{tag}_{seqtype}_contig.fasta` — FASTA of filtered contig sequences
* `{tag}_{seqtype}_clonotypes.csv` — Clonotype table with columns `clonotype_id`, `cdr3s_aa`, `frequency`, and `proportion`
* `{tag}_metrics.json` — VDJ annotation metrics (e.g. `Cells With Productive V-J Spanning Pair`, `Cells With TRA Contig`, etc.)