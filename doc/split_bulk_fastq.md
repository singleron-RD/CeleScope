# split_bulk_fastq

`celescope utils split_bulk_fastq` splits paired-end bulk FASTQ files (R1 and R2) according to the well barcodes defined by the `--well_sample` file. It generates separate `{sample}_{barcode}_R1.fastq.gz` and `{sample}_{barcode}_R2.fastq.gz` files for each sample.

This tool supports both the `bulk_rna` and `bulk_vdj` workflows. The barcode pattern is determined by the selected `--assay` and `--chemistry` parameters.

## Usage

```bash
celescope utils split_bulk_fastq \
  --assay {bulk_rna|bulk_vdj} \
  --well_sample well_sample.tsv \
  --fq1 R1.fq.gz \
  --fq2 R2.fq.gz \
  --outdir ./split_out
```

**`--assay`**  
Required. The bulk assay type. Choose `bulk_rna` or `bulk_vdj`. This determines how the barcode is parsed from R1.

**`--well_sample`**  
Required. Path to the tab-separated file mapping well numbers to sample names.

**`--fq1`**  
Required. R1 FASTQ file path. Multiple files are separated by commas.

**`--fq2`**  
Required. R2 FASTQ file path. Multiple files are separated by commas.

**`--outdir`**  
Required. Output directory for the split FASTQ files.

**`--chemistry`**  
Chemistry name. Default is `auto`. For available chemistries, see [chemistry.md](./chemistry.md).

**`--pattern`**  
Custom pattern of R1 reads, used when `--chemistry customized` is selected. For example: `C8L16C8L16C8L1U12T18`.

**`--whitelist`**  
Custom cell barcode whitelist file path, used when `--chemistry customized` is selected. One cell barcode per line.

## Output

For each sample listed in `--well_sample`, the following files are written to `--outdir`:

- `{sample}_barcode_R1.fastq.gz`: R1 reads assigned to this sample.
- `{sample}_barcode_R2.fastq.gz`: R2 reads assigned to this sample.

Reads with barcodes that are not present in the whitelist or not included in `--well_sample` are discarded.

