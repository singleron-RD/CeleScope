## Introduction

`multi_pathseq` is a CeleScope pipeline for single-cell microbial detection. It aligns sequencing reads that fail to map to the host genome against a microbial reference database using [GATK PathSeq](https://gatk.broadinstitute.org/hc/en-us/articles/360035889911--How-to-Run-the-Pathseq-pipeline), then applies the [INVADEseq](https://pmc.ncbi.nlm.nih.gov/articles/PMC10790651/) filtering strategy to produce genus-level UMI count matrices that can be merged with scRNA-Seq expression data for integrative analysis.

Pipeline steps: `sample → starsolo → pathseq → count_pathseq → analysis_pathseq`

---

## Prerequisites

1. **Download GATK PathSeq reference files** 

Users can download pre-built host and microbe reference files from Google Cloud bucket using [gsutil](https://cloud.google.com/storage/docs/gsutil_install#linux).
```bash
gsutil -m cp -r "gs://gatk-best-practices/pathseq/resources" .
```
2. **Complete scRNA-Seq analysis** with `multi_rna`. The pathseq pipeline depends on the cell barcodes and t-SNE results (accessed via `match_dir`, 4th column of the mapfile).

---

## Usage

Create a shell script `run.sh`:

```bash
multi_pathseq \
  --mapfile ./pathseq.mapfile \
  --genomeDir ./mmu \
  --filter_bwa_image ./pathseq/resources/mm10.fasta.img \
  --kmer_file ./pathseq/resources/mm10.bfi \
  --microbe_bwa_image ./pathseq/resources/pathseq_microbe.fa.img \
  --microbe_dict ./pathseq/resources/pathseq_microbe.dict \
  --microbe_taxonomy_file ./pathseq/resources/pathseq_taxonomy.db \
  --mod shell \
  --thread 8
```

Run the pipeline:

```bash
sh run.sh
bash ./shell/{sample}.sh
```

> **Note:** `./shell/{sample}.sh` must be executed from the working directory, not inside the `shell/` subdirectory.

### Arguments

`--mapfile`
Required. Tab-delimited file with 4 columns per line: (1) FASTQ prefix, (2) FASTQ directory, (3) sample name (output prefix), (4) `match_dir` from a completed `multi_rna` run.
```bash
$ cat ./pathseq.mapfile
fastq_prefix1	fastq_dir1	sample1	/path/to/rna/sample1
fastq_prefix2	fastq_dir2	sample1	/path/to/rna/sample1
fastq_prefix3	fastq_dir1	sample2	/path/to/rna/sample2
```

`--genomeDir`
Required. Host genome directory from `celescope rna mkref`.

`--filter_bwa_image`
Required. Host genome BWA image file (e.g. `mm10.fasta.img`). Included in GATK pre-built resources.

`--kmer_file`
Required. Host genome kmer file (e.g. `mm10.bfi`). Included in GATK pre-built resources.

`--microbe_bwa_image`
Required. Microbial reference BWA image file. Included in GATK pre-built resources.

`--microbe_dict`
Required. Microbial reference sequence dictionary. Included in GATK pre-built resources.

`--microbe_taxonomy_file`
Required. NCBI taxonomy database (`pathseq_taxonomy.db`). Included in GATK pre-built resources.

`--min_score_identity`
Default: `0.7`. Minimum read-to-reference identity for PathSeq alignment (`--min-score-identity`).

`--min_clipped_read_length`
Default: `60`. Minimum read length after soft-clipping for PathSeq alignment (`--min-clipped-read-length`).

`--mod`
Default: `sjm`. Script type: `shell` or [Simple Job Manager](https://github.com/StanfordBioinformatics/SJM).

`--thread`
Default: `16`. Number of threads.

---

## Pipeline Steps

`sample` and `starsolo` are shared CeleScope steps. For `starsolo` arguments and output, see [multi_rna.md](./multi_rna.md). The pathseq pipeline runs starsolo with `--outSAMunmapped Within` appended to preserve non-host reads in the BAM.

### pathseq

Extract unmapped reads from the STARsolo BAM file, downsample to a maximum of 10 million reads using Picard DownsampleSam (following the default setting in the GATK PathSeq [WDL workflow](https://github.com/broadinstitute/gatk/blob/master/scripts/pathseq/wdl/pathseq_pipeline.wdl#L187) to limit computational resource requirements), add read groups using Picard, and run gatk PathSeqPipelineSpark against the microbial reference database.

Outputs:
- `{sample}_pathseq.bam` — Reads successfully aligned to microbial references. Each carries a `YP` tag listing the NCBI taxIDs of aligned species.
- `{sample}_pathseq_score.txt` — GATK PathSeq taxonomy score summary (taxIDs, lineage names, taxonomic level, scores).
- `{sample}_unmapped.bam` — Reads that did not map to the host genome.

### count_pathseq

Counts genus-level UMI per cell using the [INVADEseq](https://pmc.ncbi.nlm.nih.gov/articles/PMC10790651/) filtering strategy:

1. Link each unmapped read to its cell barcode (`CB`) and UMI (`UB`).
2. Match reads between `{sample}_unmapped.bam` and `{sample}_pathseq.bam` (must have both `YP` and `AS` tags).
3. For each `(CB, UMI)`, keep only the alignment with the highest `AS` score.
4. Resolve all `YP` taxIDs of each retained UMI to genus names via `pathseq_score.txt`. If they all resolve to the same genus, the UMI is kept; otherwise discarded.
5. Aggregate unambiguous UMIs by genus × cell barcode to build two matrices.

Metrics (added to HTML report):

| Metric | Description |
|---|---|
| Number of Reads in Cells | Unmapped reads carrying valid `CB` and `UB`. |
| Number of Cell Reads mapping to Microbe | Reads from cells with a `YP` tag in the PathSeq BAM. |
| Number of Cell UMIs Mapping to Microbe | Unique (CB, UMI) pairs among microbe-mapping reads. |
| Number of Unambiguous UMIs | UMIs retained after the genus-consistency filter. |
| Estimated Saturation | `1 - n_unambiguous_UMIs / n_cell_microbe_reads`. |

### analysis_pathseq

Merge filtered UMI matrix with scRNA-Seq t-SNE coordinates from `match_dir`. Generates in the HTML report:
- A genus summary table (total UMI per genus, number of positive cells).
- t-SNE plot colored by scRNA-Seq Leiden cluster.
- t-SNE plot with dropdown for the top 10 genera (sorted by total UMI), colored by `log2(UMI count + 1)`.

---

## Output Files

Both matrices have **rows = genus names, columns = cell barcodes**:

- `outs/{sample}_raw_UMI_matrix.tsv.gz` — All barcodes (background + cells) with detected microbial reads. For diagnostic use.
- `outs/{sample}_filtered_UMI_matrix.tsv.gz` — Cell-matched matrix aligned to scRNA-Seq filtered cell barcodes. Cells with no microbial expression are retained as 0 to maintain column alignment. Primary matrix for downstream analysis.

An example integrative analysis with `Seurat` is available at [GitHub](https://github.com/singleron-RD/analysis_guide/blob/main/PathSeq_analysis.ipynb) or [Gitee](https://gitee.com/singleron-rd/analysis_guide/blob/main/PathSeq_analysis.ipynb).

---

## Reference

- [INVADEseq, Nat Protoc.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10790651/)
- [GATK PathSeq Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
- [GATK PathSeq Pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035889911--How-to-Run-the-Pathseq-pipeline)