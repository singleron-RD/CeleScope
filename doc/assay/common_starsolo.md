## starsolo

`starsolo` is a shared CeleScope step built on [STARsolo](https://github.com/alexdobin/STAR). It aligns R2 reads to a `mkref` genome, extracts the cell (or well) barcode and UMI from R1 according to the chemistry pattern, and collapses reads into a gene × barcode UMI count matrix. Every assay that uses it (`multi_rna`, `multi_bulk_rna`, `multi_citeseq`, ...) calls it with the same arguments; the assay-specific behaviour is described in each assay's document.

The step performs one functional run: it resolves chemistry and pattern, aligns and counts with STARsolo, post-processes the matrices, and summarizes the results into the HTML report.

### Process

#### 1. Chemistry and pattern resolution

The sample chemistry is resolved (see [common_sample.md](./common_sample.md)), and its R1 pattern is split into segments and matched with the corresponding barcode whitelist. The pattern letters are:

* `C` — cell (or well) barcode
* `L` — linker (common sequence)
* `U` — UMI
* `T` — poly-T

Depending on how many barcode segments the pattern has, the barcode and UMI positions are passed to STARsolo in one of three ways:

* **A single `C` segment**: the barcode is at a fixed offset, and the positions are given as separate start/length values. In this mode barcode matching is always strict (1 mismatch).
* **Multiple `C` segments**: each segment has its own position range, allowing linkers between barcode parts.
* **`GEXSCOPE-V3` and `flv_rna-V2`**: the barcode parts are not at a fixed offset in these libraries, so the barcode is described as three adjacent segments and the two V3 linker sequences are supplied as an adapter with one allowed mismatch. This is a workaround for [STAR issue #1607](https://github.com/alexdobin/STAR/issues/1607#issuecomment-1210096785).

The whitelist is handed to STARsolo whenever the chemistry provides one; it can be a plain text file or a gzip-compressed file, in which case a decompression stream is used, and no whitelist at all is passed when the chemistry has none.

#### 2. Alignment and counting

R2 is the genomic read and R1 the barcode read, matching STAR's convention; gzipped FASTQ files are decompressed on the fly. The following defaults apply to every assay:

* per-barcode read statistics are enabled, and the barcode read length is not checked;
* the 3′ adapter from `--adapter_3p` is clipped, and `--outFilterMatchNmin` controls the minimum matched bases;
* SAM attributes are written only when a BAM is requested;
* `--STAR_param` is applied last, so any parameter given there overrides the default of the same key.

#### 3. Feature handling

`--soloFeatures` controls which counting features STARsolo produces. The list is reordered so that `--report_soloFeature` comes first, because STARsolo uses the first feature for the BAM; `GeneFull_Ex50pAS` is always included because the intronic metrics depend on it; and `Velocyto` is added automatically when `--create_loom` is set.

#### 4. Matrix post-processing

After STAR finishes, the raw and filtered matrices are gzip-compressed in place (already-compressed files are skipped), and their directories are made world-readable. The step stops with an error if the filtered matrix is missing, the most common cause being a truncated FASTQ file.

#### 5. Q30 of barcode and UMI

The barcode and UMI quality is measured on the first 1,000,000 R1 reads of every R1 file: barcode bases and UMI bases are collected separately, and the fraction of bases with Phred score ≥ 30 is reported. If more than 10,000 reads are available, only reads 10,001–1,000,000 are counted; otherwise the first 10,000 reads are used.

#### 6. Loom output

With `--create_loom`, a `.loom` file is written for velocity analysis: the spliced, unspliced and ambiguous counts plus the gene counts are combined, restricted to the cells of the filtered matrix, and each cell ID is prefixed with the sample name.

#### 7. Report metrics

The report receives three groups of metrics:

* **Alignment metrics** — from the per-barcode read statistics of STARsolo. Region fractions are computed over valid reads or over mapped reads, as appropriate.
* **Cell metrics and barcode rank plot** — the number of cells, fraction of reads in cells, mean used reads per cell, median UMI and median genes per cell, and a log-log plot of barcodes ranked by UMI count in which called cells and background barcodes are highlighted separately.
* **Demultiplexing metrics** — raw reads, valid read fraction, corrected barcodes and Q30 values, reported as run-level numbers.

A per-barcode table (`counts.tsv`) is also written, sorted by UMI count, giving each barcode's UMI, the number of counted UMI reads and whether it was called as a cell. It is used to mark the called cells in the barcode rank plot.

### Arguments

`--chemistry`
Default: `auto`. Chemistry used for pattern and whitelist resolution. See [chemistry.md](../chemistry.md) for the full list.

`--pattern`
Optional. R1 read pattern, e.g. `C8L16C8L16C8L1U12T18`. The number after each letter is the number of bases; see [Chemistry and pattern resolution](#1-chemistry-and-pattern-resolution) for the meaning of `C`, `L`, `U` and `T`.

`--whitelist`
Optional. Cell barcode whitelist file(s), separated by whitespace, one barcode per line.

`--adapter_3p`
Default: `AAAAAAAAAAAA`. Adapter sequence(s) clipped from the 3′ end. Multiple sequences are separated by space.

`--genomeDir`
Required. Genome directory generated by `mkref` (see [multi_rna.md](multi_rna.md)).

`--outFilterMatchNmin`
Default: `50`. Alignment is output only if the number of matched bases is at least this value.

`--soloCellFilter`
Default: `EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000`. Same as the argument in STARsolo. Assays that do not perform cell calling set this to `None` (see each assay's document); their `filtered` matrix then simply contains all barcodes with at least one UMI.

`--limitBAMsortRAM`
Default: `32000000000`. Same as the argument in STARsolo.

`--STAR_param`
Additional STAR parameters, applied last so they override the default of the same key.

`--outSAMtype`
Default: `BAM SortedByCoordinate`. Same as the argument in STARsolo. Set to `None` to skip the BAM; downstream steps that read the BAM then cannot run.

`--SAM_attributes`
Additional attributes appended to the default SAM attribute list `NH HI nM AS CR UR CB UB GX GN`.

`--soloFeatures`
Default: `GeneFull_Ex50pAS Gene`. Same as the argument in STARsolo.

`--soloCBmatchWLtype`
Default: `EditDist_2`. Same as the argument in STARsolo. Note that `EditDist_2` only works when the barcode is described as multiple segments; the single-segment mode always uses strict matching.

`--report_soloFeature`
Default: `GeneFull_Ex50pAS`. Which feature of `soloFeatures` is used for the HTML report, the `outs` matrix directory and the BAM.

`--create_loom`
If set, create a `.loom` file, which requires `Velocyto` in `soloFeatures` (added automatically).

`--thread`
Default: `16`. Threads passed to STAR.

`--fq1` / `--fq2`
Required. R1 and R2 FASTQ files; multiple files of the same library are separated by commas. The number of R1 and R2 files must match, otherwise the step exits with an error.

### Output

`{outdir}/{sample}_Aligned.sortedByCoord.out.bam` — coordinate-sorted BAM with `CB`, `UB`, `GX`, `GN` and `NH` tags, used by `count`, `cells` and any BAM-splitting step.

`{outdir}/{sample}_Solo.out/{report_soloFeature}/raw/` — raw matrix in 10x Genomics format, gzip-compressed: `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`.

`{outdir}/{sample}_Solo.out/{report_soloFeature}/filtered/` — filtered (cell-called) matrix, same format and compression.

`{outdir}/{sample}_Solo.out/{report_soloFeature}/Summary.csv` — STARsolo summary statistics, including sequencing saturation and Q30 of RNA reads.

`{outdir}/{sample}_Solo.out/{report_soloFeature}/CellReads.stats` — per-barcode read statistics from which the alignment metrics are derived.

`{outdir}/{sample}_Solo.out/{report_soloFeature}/counts.tsv` — per-barcode `UMI`, `countedU` and `mark` (`CB` / `UB`) table sorted by UMI.

`{outdir}/{sample}.loom` — only with `--create_loom`.

All files under `{sample}_Solo.out/` are made world-readable.

### Metrics (added to HTML report)

Alignment and barcode metrics, derived from `CellReads.stats`:

| Metric Name | Description |
|---|---|
| Genome | Genome name read from the `--genomeDir` configuration |
| Reads Mapped to Unique Loci | Reads mapping uniquely to the genome, over valid reads |
| Reads Mapped to Multiple Loci | Reads mapping to multiple loci, over valid reads |
| Reads Mapped Uniquely to Transcriptome | Reads counted by STARsolo (one gene per read); these reads are used for UMI counting |
| Mapped Reads Assigned to Exonic Regions | Fraction of mapped reads overlapping exons |
| Mapped Reads Assigned to Intronic Regions | Fraction of mapped reads overlapping introns |
| Mapped Reads Assigned to Intergenic Regions | Fraction of mapped reads that cannot be assigned to a gene |
| Mapped Reads Antisense to Gene | Fraction of mapped reads on the opposite strand of a gene |

Cell-level metrics, from `Summary.csv` and the filtered matrix:

| Metric Name | Description |
|---|---|
| Estimated Number of Cells | Number of barcodes considered cell-associated, i.e. present in the filtered matrix |
| Fraction Reads in Cells | `Fraction of Unique Reads in Cells` from `Summary.csv`: uniquely-mapped-to-transcriptome reads with cell-associated barcodes |
| Mean Reads per Cell | Valid reads divided by the estimated number of cells |
| Mean Used Reads per Cell | Mean number of uniquely-mapped-to-transcriptome reads per cell-associated barcode |
| Median UMI per Cell | Median UMI count per cell-associated barcode |
| Total Genes | Number of genes with at least one UMI count in any cell |
| Median Genes per Cell | Median number of detected genes per cell-associated barcode |
| Saturation | Sequencing saturation reported by STARsolo |
| Q30 of RNA Reads | Fraction of RNA read bases with quality score ≥ 30, over all reads |

Demultiplexing metrics:

| Metric Name | Description |
|---|---|
| Raw Reads | Total reads from the FASTQ files |
| Valid Reads | Fraction of reads with a valid barcode and UMI |
| Corrected Barcodes | Fraction of valid reads whose barcode was corrected to a whitelist barcode 1 Hamming distance away |
| Q30 of Barcodes | Fraction of barcode bases with quality score ≥ 30, from reads 10,001–1,000,000 |
| Q30 of UMI | Fraction of UMI bases with quality score ≥ 30, from reads 10,001–1,000,000 |
| Q30 of RNA Reads | Fraction of RNA read bases with quality score ≥ 30, over all reads |

`--report_soloFeature` and `--soloCellFilter` (hidden) are also recorded in the report.

### HTML Report Plots

| Plot | Description |
|---|---|
| Barcode rank plot | Log-log plot of barcodes ranked by UMI count, split into called cells and background by the `CB` / `UB` mark |

### Notes

* The BAM is produced with the first element of `soloFeatures`, which is why `--report_soloFeature` is moved to the front of the list.
* `GeneFull_Ex50pAS` is always part of `soloFeatures`, even if the user list omits it, because the intronic and antisense metrics require it.
* `--outSAMtype None` disables the BAM, the SAM attribute output and every downstream step that reads the BAM.