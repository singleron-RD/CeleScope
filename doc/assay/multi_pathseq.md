## **Usage**  

### Obtain pre-built reference files from GATK
Users can download pre-built host and microbe reference files from Google Cloud bucket.
[install gsutil](https://cloud.google.com/storage/docs/gsutil_install#linux)
```
gsutil -m cp -r "gs://gatk-best-practices/pathseq/resources" .
```

**Reference**
1.https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle

2.https://gatk.broadinstitute.org/hc/en-us/articles/360035889911--How-to-Run-the-Pathseq-pipeline

### Generate scripts for each sample
In your working directory, create a shell script named `run.sh` with the following content:  

```bash
multi_pathseq \
 --mapfile mapfile \
 --genomeDir ./mmu \
 --filter_bwa_image ./mm10.fasta.img \
 --kmer_file ./mm10.bfi \
 --microbe_bwa_image ./pathseq_microbe.fa.img \
 --microbe_dict ./pathseq_microbe.dict \
 --microbe_taxonomy_file ./pathseq_taxonomy.db \
 --mod shell \
 --thread 8
```

#### **Parameter Descriptions:**  
- `--mapfile` (**Required**): A tab-delimited text file with at least three columns, where each row represents a paired-end FASTQ file entry.  

  **Column structure:**  
  1st column: FASTQ file prefix  
  2nd column: Directory path containing the FASTQ files  
  3rd column: sample name(used as the prefix for output matrix file)  
  4th column: The single cell rna directory after running CeleScope

  **Example:**  
  Sample1 has two paired-end FASTQ files located in different directories (`fastq_dir1` and `fastq_dir2`), while Sample2 has one paired-end FASTQ file in `fastq_dir1`.  

  ```bash
  $ cat ./my.mapfile
  fastq_prefix1  fastq_dir1  sample1    rna/sample1
  fastq_prefix2  fastq_dir2  sample1    rna/sample1
  fastq_prefix3  fastq_dir1  sample2    rna/sample2

  $ ls fastq_dir1
  fastq_prefix1_1.fq.gz  fastq_prefix1_2.fq.gz
  fastq_prefix3_1.fq.gz  fastq_prefix3_2.fq.gz

  $ ls fastq_dir2
  fastq_prefix2_1.fq.gz  fastq_prefix2_2.fq.gz
  ```  

- `--genomeDir` (**Required**): Path to the genome directory created using `celescope rna mkref`.  
- `--filter_bwa_image`, `--kmer_file`, `--microbe_bwa_image`, `--microbe_dict`, `--microbe_taxonomy_file` (**Required**): These files are included in the pre-built reference files.

- `--mod`: Specifies the script format:  
  - `shell` bash script
  - `sjm` [Simple Job Manager](https://github.com/StanfordBioinformatics/SJM) script
- `--thread`: Number of threads to use (maximum: 20).  

### Run the analysis  
After executing `sh run.sh`, a `shell/` directory will be created, containing `{sample}.sh` scripts.  

Start the analysis by running:  
```bash
bash ./shell/{sample}.sh
```  
> **Note:** The `{sample}.sh` script must be executed from the working directory, not from inside the `shell/` directory.  

---

## Main Output Files

- `outs/{sample}_raw_UMI_matrix.tsv.gz`  
  - UMI matrix, where rows represent microbes and columns correspond to cell barcodes.


