## `celescope rna cells`: additional matrix filtering

After completing the `multi_rna` workflow, you can run `celescope rna cells` on the existing results to re-perform cell calling. This step **does not require the BAM file** — it operates only on matrix-level artifacts from the quantification step. You can safely archive BAMs before running it.

### What it reads

| File | Purpose |
|---|---|
| `{sample}/outs/raw/` | Unfiltered raw UMI matrix (all whitelist barcodes) |
| `{sample}/outs/counts.tsv` | Per-barcode summary table (UMI, read counts, etc.) |
| Pipeline slot files | Internal metrics (e.g. Saturation) from the original run |

> **Important**: Because it reads pipeline slot files left by the original `multi_rna` workflow, `celescope rna cells` must run from the **same working directory** where `multi_rna` was executed, or you must point to that directory with `--root_dir`. You cannot copy just a bare `outs/` directory to another location and run it in isolation.


### Filtering modes (mutually exclusive)

The three cell-selection arguments are mutually exclusive. Priority order (highest first):

1. `--barcode` — use a user-provided list of barcodes (one per line)
2. `--force_cells` — take the top N barcodes ranked by UMI count
3. `--soloCellFilter` — run STARsolo's EmptyDrops-like filtering

### Additional QC filters

These can be combined with any of the three modes above, and are applied **after** cell selection:

| Argument | Default | Description |
|---|---|---|
| `--min_gene` | `0` (disabled) | Drop cells with fewer than N detected genes |
| `--max_mito` | `1.0` (disabled) | Drop cells whose mitochondrial UMI fraction exceeds this value. **Requires `--genomeDir`** |

### Output behavior

- On the **first run**, the original `outs/filtered/` matrix and `{sample}_report.html` are renamed with a `default_` prefix (e.g. `default_filtered/`, `default_{sample}_report.html`).
- On **subsequent runs**, the `default_` prefix is not touched — the step is idempotent. You can safely rerun with different filtering parameters and never lose the original results.
- A new `outs/filtered/` matrix and `{sample}_report.html` are generated each time.

### Arguments

| Argument | Required | Description |
|---|---|---|
| `--sample` | ✅ | Sample name. The step looks under `{sample}/outs/` for existing results. |
| `--root_dir` | No (default `./`) | Absolute path to the `multi_rna` working directory. Set this if submitting from a job scheduler that launches from `/tmp` or another directory. |
| `--outdir` | ✅ (**but ignored**) | Required by the argparse template but internally overwritten to `{sample}/cells`. Any placeholder value works. |
| `--force_cells` | No | Force this many cells (top N by UMI). |
| `--barcode` | No | Path to a user-defined barcode file (one per line). |
| `--soloCellFilter` | No | STARsolo EmptyDrops-like filter string. |
| `--min_gene` | No | Minimum detected genes per cell. |
| `--max_mito` | No | Maximum mitochondrial UMI fraction per cell. Requires `--genomeDir`. |
| `--genomeDir` | Conditional | Required only when `--max_mito < 1.0`. |

### Example 1: Force cell number with QC filters (run from the original working directory)

```bash
# cd to the directory where multi_rna was originally executed
cd /path/to/your/working_dir

celescope rna cells \
  --force_cells 10000 \
  --min_gene 200 \
  --max_mito 0.05 \
  --genomeDir /path/to/mmu_or_hs_genome \
  --outdir dummy \
  --sample test1
```

### Example 2: Force cell number from a different directory using `--root_dir`

```bash
celescope rna cells \
  --force_cells 10000 \
  --root_dir /path/to/your/working_dir \
  --outdir dummy \
  --sample test1
```

This is useful when submitting through a job scheduler that starts from `/tmp` or another scratch directory. Internally, the step does `os.chdir(root_dir)`, so you only need to point to the absolute path of the original working directory.

### Example 3: Adjust cell-calling parameters with `--soloCellFilter`

For more details about STARsolo EmptyDrops-like filtering, please refer to:
[https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#emptydrop-like-filtering](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#emptydrop-like-filtering)

```bash
# In CeleScope v1.*, the default FDR for cell calling is 0.01.
# In CeleScope v2.*, the default FDR is 0.001.
# The following parameters set the FDR to 0.01:

celescope rna cells \
  --soloCellFilter "EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000" \
  --root_dir /path/to/your/working_dir \
  --outdir dummy \
  --sample test1
```