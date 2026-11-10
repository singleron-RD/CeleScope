## cutadapt

`cutadapt` is the adapter-trimming step shared by most CeleScope assays ( `multi_capture_virus`,  `multi_snp`, `multi_citeseq`, etc.). It consumes the demultiplexed R2 FASTQ written by `barcode` (see [common_barcode.md](./common_barcode.md)) and produces the trimmed R2 FASTQ used by all downstream steps. The step is implemented in `celescope/tools/cutadapt.py` and is invoked with the same arguments regardless of the assay.

Trim sequencing adapters from R2 reads with [cutadapt](https://cutadapt.readthedocs.io/).

### Process

Run cutadapt on the demultiplexed R2 reads to remove the default adapter (`-a polyA=A{18}`, i.e. 18 A bases). `--nextseq-trim` (default `20`) is always passed for two-color chemistry quality trimming: NextSeq / NovaSeq encode a "dark cycle" (no detected color) as a `G`, and such dark cycles occur when sequencing falls off the end of the fragment, producing a run of high-quality but incorrect `G` calls at the 3' end.

Because cutadapt allows partial matches between read and adapter, short matches can occur by chance and trim bases erroneously. The `--overlap` parameter requires that at least `overlap` bases match between adapter and read, reducing false trimming.

When `--skip_cutadapt` is set, the input FASTQ is copied to the output FASTQ unchanged.

### Arguments

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

### Output

`{outdir}/{sample}_cutadapt.json` — cutadapt JSON log.

`{outdir}/{sample}_clean_2.fq` — adapter-trimmed R2 FASTQ. This file is consumed by the downstream steps (for example `star_virus` when `--not_consensus` is used, or `consensus`).

### Metrics (added to HTML report)

| Metric Name | Description |
|---|---|
| Reads with Adapters | Reads with poly-A tails (and user-provided adapters, if any) that were trimmed |
| Reads too Short | Reads whose length is less than `--minimum_length` after trimming |
| Base Pairs Quality-Trimmed | Base pairs removed from the read end whose quality is below the `--nextseq-trim` threshold |