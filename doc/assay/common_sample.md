## sample

`sample` is the first step of every `multi_*` assay. It is an internal bookkeeping and sanity-check step: it does not transform reads for analysis, but it detects the single-cell chemistry, writes the metrics header of the HTML report, and produces a small debug bundle that helps diagnose barcode/linker problems before the real pipeline runs. The step is implemented in `celescope/tools/sample.py` and is invoked with the same arguments regardless of the assay.

It is step `00.sample` in each sample's output directory, ahead of `01.barcode`.

Its input is the FASTQ list assembled from the mapfile: for each sample name (3rd column), all rows sharing that name contribute one R1/R2 pair, and the pairs are passed as comma-separated lists.

### Process

1. **Chemistry detection.** `get_chemistry` resolves the chemistry for the sample. When `--chemistry auto` is given, up to the first 10,000 R1 reads are scanned and each read is matched against every chemistry's barcode whitelist (and, for `GEXSCOPE-V3` / `flv_rna-V2`, the linkers). The chemistry with the most matching reads wins. Detection fails if the winning chemistry covers less than 10% of the scanned reads; a warning is printed if it covers less than 50%. All FASTQ files of one sample must resolve to the same chemistry, otherwise the step exits with an error.

2. **Report header.** Four metrics identify the run: sample ID, assay, detected chemistry and software version.

3. **Debug subsampling.** For every input FASTQ, `--debug_skip_read` reads are skipped and then `--debug_use_read` reads are written to `{fq_stem}.subsample.gz`. Skipping is needed because the first hundred reads of an R1 FASTQ typically contain `N` bases at the beginning.

4. **Invalid-read rendering.** For the first R1 file, the same read window is re-parsed with the resolved chemistry, and every read whose barcode cannot be validated is written to an HTML page. Each observed barcode segment (`bc1`, `bc2`, `bc3`) and linker (`linker1`, `linker2`, `linker3`) is highlighted with its own background colour, with the colour legend shown at the top, so that the exact base positions responsible for a failed match can be seen at a glance.

5. **Packing.** The subsampled FASTQs and the invalid-read HTML pages are collected into a single tar archive.

### Arguments

`--chemistry`
Default: `auto`. Single-cell chemistry used for detection. Choices are the keys of `chemistry_dict`. If set explicitly, detection is skipped and the given chemistry is used. See [chemistry.md](../chemistry.md) for the list of chemistries and their barcode structures.

`--pattern`
Optional. R1 read pattern, e.g. `C8L16C8L16C8L1U12T18`, used only together with `--chemistry customized`. The number after each letter is the number of bases:
* `C` — cell barcode
* `L` — linker (common sequence)
* `U` — UMI
* `T` — poly-T

`--whitelist`
Optional. Cell barcode whitelist file(s), separated by whitespace, one barcode per line. Used together with `--chemistry customized`, or to override the chemistry's built-in whitelist.

`--linker`
Optional. Linker sequence file(s), separated by whitespace. Used together with `--chemistry customized`.

`--debug_use_read`
Default: `1000`. Number of reads rendered into the debug outputs.

`--debug_skip_read`
Default: `10000`. Number of reads skipped before the debug window, since the R1 FASTQ usually contains `N` bases in the first hundred reads.

### Output

`{outdir}/{fq_basename}.subsample.gz` — one per input FASTQ; a small subsample taken from the debug window. The name is built from the FASTQ file name with its final suffix removed, so `xxx_1.fq.gz` yields `xxx_1.fq.subsample.gz`.

`{outdir}/{fq1_basename}.invalid.html` — colour-coded rendering of the invalid-barcode R1 reads, with a colour legend for each barcode segment and linker. Only the first R1 file is rendered.

`{outdir}/{sample}_debug.tar` — archive containing the subsampled FASTQs and the invalid-read HTML pages.

No FASTQ consumed by `barcode` is produced here; `barcode` reads the original files directly.

### Metrics (added to HTML report)

| Metric Name | Description |
|---|---|
| Sample ID | Sample name, taken from the 3rd mapfile column |
| Assay | Assay name |
| Chemistry | Detected (or given) chemistry; help text links to the chemistry documentation |
| Software Version | CeleScope version |

No plots are produced by this step.