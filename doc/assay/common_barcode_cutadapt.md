## Common steps: barcode and cutadapt

`barcode` and `cutadapt` are the first two data-processing steps shared by most CeleScope assays (`multi_rna`, `multi_capture_virus`, `multi_tag`, `multi_vdj`, `multi_snp`, `multi_citeseq`, etc.). They start from the raw paired-end FASTQ files and produce a demultiplexed, adapter-trimmed R2 FASTQ that all downstream steps consume.

`barcode → cutadapt`

Both steps are implemented in `celescope/tools/` and are invoked with the same arguments regardless of the assay.

---

### barcode

Demultiplex R1 reads: extract cell barcodes and UMIs from R1, correct barcodes against a whitelist, and write corrected R2 with demultiplexed read names.

#### Process

Detect single-cell chemistry (default `auto`) and parse the R1 pattern into cell barcode (`C`), linker (`L`), UMI (`U`) and poly-T (`T`) segments against the corresponding barcode whitelist. A precomputed lookup table allows up to 1 Hamming distance from each whitelist barcode, so that observed barcode sequences can be corrected in constant time. For each paired R1–R2 read, the cell barcode and UMI are extracted from R1, the barcode is corrected through the lookup table, and the read is written to a demultiplexed R2 FASTQ whose header is reformatted as `{corrected_barcode}:{umi}:{read_index}`. Reads whose barcode cannot be corrected are discarded.

For `GEXSCOPE-V3` and `flv_rna-V2`, a per-read offset is detected first and the sequence is shifted accordingly before parsing. For `flv` chemistries the barcode is taken from the reverse complement of the reversed C segments.

#### Arguments

`--chemistry`
Default: `auto`. Single-cell chemistry. Choices: keys of `chemistry_dict` (e.g. `auto`, `GEXSCOPE-V1`, `GEXSCOPE-V2`, `GEXSCOPE-V3`, `flv`, `flv-V2`, ...). Use `--chemistry customized` together with `--pattern` and, optionally, `--whitelist` for user-defined combinations. For legacy chemistries without a built-in whitelist, the chemistry must be given explicitly. See [chemistry.md](../chemistry.md) for the full list.

`--pattern`
Optional. R1 read pattern, e.g. `C8L16C8L16C8L1U12T18`. The number after each letter is the number of bases:
* `C` — cell barcode
* `L` — linker (common sequence)
* `U` — UMI
* `T` — poly-T

`--whitelist`
Optional. Cell barcode whitelist file, one barcode per line. Needed when overriding the chemistry's built-in whitelist or when using `--chemistry customized`.

#### Output

`{outdir}/{sample}_2.fq` — demultiplexed R2 FASTQ (no `.gz` suffix). Only reads with valid or correctable cell barcodes are kept, and the header is `{corrected_barcode}:{umi}:{read_index}`. This file is the input of `cutadapt`.

#### Metrics (added to HTML report)

| Metric Name | Description |
|---|---|
| Raw Reads | Total reads from the FASTQ files |
| Valid Reads | Reads with a correct cell barcode |
| Corrected Reads | Reads whose barcode was not in the whitelist but was within one Hamming distance of it, and was corrected |
| Q30 of Barcode | Fraction of cell barcode bases with Phred quality score ≥ 30 |
| Q30 of UMI | Fraction of UMI bases with Phred quality score ≥ 30 |

---

### cutadapt

Trim sequencing adapters from R2 reads with [cutadapt](https://cutadapt.readthedocs.io/).

#### Process

Run cutadapt on the demultiplexed R2 reads to remove the default adapter (`-a polyA=A{18}`, i.e. 18 A bases). `--nextseq-trim` (default `20`) is always passed for two-color chemistry quality trimming: NextSeq / NovaSeq encode a "dark cycle" (no detected color) as a `G`, and such dark cycles occur when sequencing falls off the end of the fragment, producing a run of high-quality but incorrect `G` calls at the 3' end.

Because cutadapt allows partial matches between read and adapter, short matches can occur by chance and trim bases erroneously. The `--overlap` parameter requires that at least `overlap` bases match between adapter and read, reducing false trimming.

When `--skip_cutadapt` is set, the input FASTQ is copied to the output FASTQ unchanged.

#### Arguments

`--minimum_length`
Default: `20`. Discard processed reads shorter than this length.

`--overlap`
Default: `10`. Minimum number of bases that must match between the adapter and the read.

`--nextseq_trim`
Default: `20`. Quality trimming threshold for two-color chemistry instruments.

`--cutadapt_param`
Other cutadapt parameters, enclosed in quotation marks, e.g. `--cutadapt_param "-a p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"` to additionally trim the Illumina p5 adapter. There is no `--adapter_fasta` argument; extra adapters must be supplied through this argument.

`--skip_cutadapt`
If set, skip cutadapt and copy the input FASTQ to the output FASTQ directly.

#### Output

`{outdir}/{sample}_cutadapt.json` — cutadapt JSON log.

`{outdir}/{sample}_clean_2.fq` — adapter-trimmed R2 FASTQ. This file is consumed by the downstream steps (for example `star_virus` when `--not_consensus` is used, or `consensus`).

#### Metrics (added to HTML report)

| Metric Name | Description |
|---|---|
| Reads with Adapters | Reads with poly-A tails (and user-provided adapters, if any) that were trimmed |
| Reads too Short | Reads whose length is less than `--minimum_length` after trimming |
| Base Pairs Quality-Trimmed | Base pairs removed from the read end whose quality is below the `--nextseq-trim` threshold |

---

## Downstream usage

`{sample}_clean_2.fq` becomes the common input of the remaining pipeline. In an assay without a consensus step it is mapped directly; in `multi_capture_virus` it is first passed to `consensus` unless `--not_consensus` is given.

Other assays that reuse these two steps may add their own post-processing; refer to the corresponding assay document.