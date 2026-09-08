## Reference

- Human
```
mkdir -p /genome/vdj/human
cd /genome/vdj/human
wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TR{A,B}{V,J}.fasta
wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBD.fasta
wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IG{H,K,L}{V,J}.fasta
wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta
celescope vdj mkref human TR
celescope vdj mkref human IG
```

- Mouse
```
mkdir -p /genome/vdj/mouse
cd /genome/vdj/mouse
wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TR{A,B}{V,J}.fasta
wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TRBD.fasta
wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IG{H,K,L}{V,J}.fasta
wget https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGHD.fasta
celescope vdj mkref mouse TR
celescope vdj mkref mouse IG
```

## Usage


### Generate scripts for each sample

In your working directory, create a shell script named `run.sh`:

```bash
multi_bulk_vdj \
    --mapfile ./vdj.mapfile \
    --ref_path /genome/vdj/human/human_TR \
    --species human \
    --type TCR \
    --well_sample well_sample.tsv \
    --thread 8 \
    --mod shell
```

#### Arguments

`--mapfile`
Required. The mapfile is a tab-delimited text file containing at least three columns. Each line represents one pair of FASTQ files.

* 1st column: FASTQ file prefix
* 2nd column: FASTQ file directory path
* 3rd column: The name of the plate.

Example:

```bash
$ cat ./my.mapfile
fastq_prefix1	fastq_dir	plate1
fastq_prefix2	fastq_dir	plate1
```

`--ref_path`
Required. The path of the reference directory after running `celescope vdj mkref`.

`--species`
Default `human`. Species name.

- When `human` or `mouse` is specified, the pipeline will use the built-in auxiliary data file (`optional_file/{species}_gl.aux`) for IgBLAST.
- For other species, you need to provide a custom auxiliary data file using `--aux_file`.

**Example workflow:**
```bash
# For human/mouse samples
multi_bulk_vdj \
    --mapfile ./vdj.mapfile \
    --ref_path /genome/vdj/human/human_TR \
    --species human \
    --type TCR \
    --well_sample well_sample.tsv \
    --mod shell

# For other species (requires custom aux_file)
multi_bulk_vdj \
    --mapfile ./vdj.mapfile \
    --ref_path /genome/vdj/custom/custom_TR \
    --species custom_species \
    --type TCR \
    --well_sample well_sample.tsv \
    --aux_file /path/to/custom_gl.aux \
    --mod shell
```

`--aux_file`
Custom auxiliary data file for IgBLAST. Required when `--species` is not `human` or `mouse`. The auxiliary data file contains information about CDR3 boundaries and framework regions needed by IgBLAST for proper V(D)J annotation.

`--type`
Required. Receptor type. `TCR` or `BCR`.

`--well_sample`
Required. A TSV file containing well numbers and sample names of wells.

* 1st column: Well numbers
* 2nd column: Corresponding sample names

**96 well number (8 × 12)**

![](../images/96-well.png)

Example:
```tsv
1	control1
2	control2
3	treatment1
4	treatment2
...
```

`--chemistry` 
Default is `auto`, which automatically detects the chemistry from the FASTQ files.

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

## Steps

### mkref

- Build index for IMGT reference sequences.
- Make sure the current directory contains all V, D, J reference sequences of TRA/TRB or IGH/IGK/IGL downloaded from the IMGT website.

### sample

- Parse sample information from the mapfile.
- Auto-detect chemistry version if `--chemistry auto` is specified.

### barcode

- Demultiplex barcodes from R1 reads.
- Filter invalid R1 reads, which includes:
    - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.
    - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.
    - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
    - Low quality reads: low sequencing quality in barcode and UMI regions.

### consensus

- Consensus all the reads of the same (barcode, UMI) combinations into one read (UMI).
- It will go through the sequence residue by residue and count up the number of each type of residue (i.e. A or G or T or C for DNA) in all sequences in the alignment.
- If the following conditions are met, the consensus sequence will be the most common residue in the alignment:
    1. The percentage of the most common residue type > `threshold` (default: 0.5);
    2. Most common residue reads >= `min_consensus_read` (default: 2);
- Otherwise, an ambiguous character (N) will be added.
- Filter out consensus sequences with more than 5 N bases.


### mapping_vdj

- Align consensus sequences to IMGT (http://www.imgt.org/) database sequences using IgBLAST.
- Parse AIRR format output and filter productive rearrangements.
- Generate clonotype annotations and diversity metrics (n_clonotypes, inverse Simpson index) for each well/sample.

---

## Main output

- `outs/annotation/` This directory contains V(D)J annotations for each well/sample.
    - Each file is named `{sample}_annotation.csv`
    - Columns: barcode, well, sample, chain, v_gene, d_gene, j_gene, productive, cdr3_nt, cdr3, raw_clonotype_id, umis, percent
    - This directory can be imported into [immunarch](https://github.com/immunomind/immunarch) for downstream analysis using the following code:
    ```r
    library(immunarch)
    immdata <- repLoad("path to outs/annotation")
    ```

- `outs/clonotypes/` This directory contains clonotype information for each well/sample.
    - Each file is named `{sample}_clonotypes.csv`
    - Aggregated based on (v_gene, d_gene, j_gene, cdr3_nt)
    - Columns: well, sample, chain, v_gene, d_gene, j_gene, cdr3_nt, cdr3, raw_clonotype_id, umis, percent

- `outs/{sample}_filtered_annotations.csv` Combined V(D)J annotations across all wells/samples.

- `outs/{sample}_clonotypes.csv` Combined clonotype information across all wells/samples.

---

## Intermediate files

### mkref

- VDJ IMGT reference with index files.

### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of the read name is `{barcode}_{UMI}_{read ID}`.

### consensus

- `{sample}_consensus.fq` Fastq file after consensus.
- `{sample}_filtered_consensus.fasta` Filtered consensus sequences with <= 5 N bases.
- `{sample}_metrics.tsv` Metrics file containing UMI statistics.

### mapping_vdj

- `{plate}_airr.tsv` The alignment result of each UMI.
    - A tab-delimited file compliant with the [AIRR Rearrangement schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)
    - Contains V, D, J gene calls, CDR3 sequences, and productivity information

- `{plate}_well_metrics.tsv` Well-level metrics including:
    - well: Well number
    - sample: Sample name
    - n_clonotypes: Total number of clonotypes
    - inverse_simpson: Inverse Simpson diversity index
    - read: Total reads
    - umi: Total UMIs
    - umi_mapped: UMIs mapped to any VDJ gene
    - umi_confident: UMIs with confident VJ mapping (see Metrics section below for detailed criteria)
    - umi_confident_{chain}: UMIs mapped to specific chains (TRA/TRB for TCR, IGH/IGK/IGL for BCR)

---

## Metrics

| Metric | Description |
|--------|-------------|
| UMIs Mapped to Any VDJ Gene | UMIs mapped to any germline VDJ gene segments |
| UMIs Mapped Confidently to VJ Gene | UMIs with productive rearrangement mapped to VJ gene pairs and with correct junction sequence. A UMI is considered "confidently mapped" when it meets ALL of the following criteria:<br>1. `productive` field is `T` (true productive rearrangement)<br>2. CDR3 nucleotide sequence (`junction`) is not empty<br>3. No ambiguous bases (`N`) in CDR3 nucleotide sequence<br>4. No stop codons (`*`) in CDR3 amino acid sequence (`junction_aa`)<br>5. CDR3 amino acid length > 5<br>6. CDR3 amino acid starts with Cysteine (`C`)<br>7. Locus is in the expected chain set (TRA/TRB for TCR, IGH/IGK/IGL for BCR) |
| UMIs Mapped Confidently to {chain} | UMIs mapped confidently to specific chains (TRA/TRB for TCR, IGH/IGK/IGL for BCR). The same confident criteria as "UMIs Mapped Confidently to VJ Gene" apply, but counted separately for each chain. |
| Fraction of Reads in Wells | Fraction of total reads assigned to wells in the well_sample file; lower value indicates wrong well_sample file or high ambient contamination |
| Filtered UMI Counts | Consensus UMIs after filtering out sequences with > 5 N bases |