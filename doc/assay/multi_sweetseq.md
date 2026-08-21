## Usage

Before running `multi_sweetseq`, you need to run scRNA-Seq data with CeleScope first to generate the matched `matched_dir`.

```
multi_sweetseq \
    --mapfile ./sweetseq.mapfile \
    --barcode_fasta celescope/data/sweetseq/sweet_tag_barcode.fasta \
    --linker_fasta celescope/data/sweetseq/sweet_tag_linker.fasta \
    --fq_pattern L23C15 \
    --mod shell
```

- `--mapfile` Mapfile is a tab-delimited text file with at least three columns. Each line of mapfile represents paired-end fastq files.

    1st column: Fastq file prefix.  
    2nd column: Fastq file directory path.  
    3rd column: Sample name, which is the prefix of all output files.  
    4th column: `matched_dir` from CeleScope scRNA-Seq analysis (required for `sweetseq`).

    Example

    Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.

    ```
    $cat ./my.mapfile
    fastq_prefix1	fastq_dir1	sample1	/path/to/matched_dir
    fastq_prefix2	fastq_dir2	sample1	/path/to/matched_dir
    fastq_prefix3	fastq_dir1	sample2	/path/to/matched_dir2

    $ls fastq_dir1
    fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
    fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

    $ls fastq_dir2
    fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
    ```

- `--barcode_fasta` Required. Tag barcode fasta file. It will check the mismatches between tag barcode sequence in R2 reads with all tag barcode sequences in barcode_fasta. It will assign read to the tag with mismatch <= threshold. If no such tag exists, the read is classified as invalid.

    The mismatch threshold is chosen automatically based on the number of tag barcodes:
    - 2 mismatches when the number of tag barcodes is <= 10000.
    - 1 mismatch when the number of tag barcodes is > 10000.

    You can find the barcode fasta file under `celescope/data/sweetseq`.

- `--linker_fasta` Optional. If provided, it will check the mismatches between linker sequence in R2 reads with all linker sequences in linker_fasta. If the minimum mismatch is not < len(linker) / 10 + 1, the read is classified as invalid.

- `--fq_pattern` R2 read pattern. The number after the letter represents the number of bases. CLindex is `L25C15` and sweetseq is `L23C15`.
    - `L` linker(common sequences)
    - `C` tag barcode.

- `--mod` Which type of script to generate, `sjm` or `shell`.

## Features

The sweetseq pipeline contains the following steps:

### sample

- Generate sample-level metadata and summary.

### barcode

- Demultiplex cell barcodes from R1 reads according to the selected chemistry or user-defined `--pattern`.
- Correct cell barcode sequences within one mismatch of the whitelist.
- Filter R1 reads that do not match the expected barcode/linker/UMI/polyT pattern.
- Output valid R2 reads with the read name formatted as `{barcode}:{UMI}:{read ID}`.

### mapping_tag

- Align R2 reads to the tag barcode fasta.
- If `--linker_fasta` is provided, validate the linker sequence in R2 reads.

### count_tag

- Assign tag to each scRNA-Seq cell barcode and summarize UMI counts.

### analysis_tag

- Combine scRNA-Seq clustering information with tag assignment.
- Generate t-SNE plots colored by tag UMI counts.

## Output files

### sample

- Sample-level log and metadata files under `00.sample/{sample}/`.

### barcode

- `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of the read name is `{barcode}:{UMI}:{read ID}`.

### mapping_tag

- `{sample}_read_count.tsv` tab-delimited text file with 4 columns.

    `barcode` cell barcode  
    `tag_name`  tag name in barcode_fasta  
    `UMI`   UMI sequence  
    `read_count` read count per UMI

- `{sample}_invalid_barcode.tsv` tab-delimited text file with 2 columns.

    `tag_barcode` tag barcodes that do not match with any sequence in `--barcode_fasta`.  
    `read_count` invalid tag barcode read counts

### count_tag

- `{sample}_umi_tag.tsv`

    `first column` cell barcode  
    `last column`  assigned tag  
    `columns between first and last` UMI count for each tag

### analysis_tag

- `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster information.
