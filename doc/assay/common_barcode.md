## barcode

`barcode` is the first data-processing step shared by most CeleScope assays (`multi_rna`, `multi_capture_virus`, `multi_tag`, `multi_vdj`, `multi_snp`, `multi_citeseq`, etc.), following the `sample` step (see [common_sample.md](./common_sample.md)) which detects the chemistry and writes the debug bundle. It starts from the raw paired-end FASTQ files and writes a demultiplexed R2 FASTQ that the downstream steps consume. The step is implemented in `celescope/tools/barcode.py` and is invoked with the same arguments regardless of the assay.

Demultiplex R1 reads: extract cell barcodes and UMIs from R1, correct barcodes against a whitelist, and write corrected R2 with demultiplexed read names.

### Process

Detect single-cell chemistry (default `auto`) and parse the R1 pattern into cell barcode (`C`), linker (`L`), UMI (`U`) and poly-T (`T`) segments against the corresponding barcode whitelist. A precomputed lookup table allows up to 1 Hamming distance from each whitelist barcode, so that observed barcode sequences can be corrected in constant time. For each paired R1–R2 read, the cell barcode and UMI are extracted from R1, the barcode is corrected through the lookup table, and the read is written to a demultiplexed R2 FASTQ whose header is reformatted as `{corrected_barcode}:{umi}:{read_index}`. Reads whose barcode cannot be corrected are discarded.

For `GEXSCOPE-V3` and `flv_rna-V2`, a per-read offset is detected first and the sequence is shifted accordingly before parsing. For `flv` chemistries the barcode is taken from the reverse complement of the reversed C segments.

### Arguments

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

### Output

`{outdir}/{sample}_2.fq` — demultiplexed R2 FASTQ (no `.gz` suffix). Only reads with valid or correctable cell barcodes are kept, and the header is `{corrected_barcode}:{umi}:{read_index}`. In assays with a trimming step this file is the input of `cutadapt` (see [common_cutadapt.md](./common_cutadapt.md)); assays without trimming consume it directly.

### Metrics (added to HTML report)

| Metric Name | Description |
|---|---|
| Raw Reads | Total reads from the FASTQ files |
| Valid Reads | Reads with a correct cell barcode |
| Corrected Reads | Reads whose barcode was not in the whitelist but was within one Hamming distance of it, and was corrected |
| Q30 of Barcode | Fraction of cell barcode bases with Phred quality score ≥ 30 |
| Q30 of UMI | Fraction of UMI bases with Phred quality score ≥ 30 |