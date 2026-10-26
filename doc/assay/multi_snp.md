## Introduction

`multi_snp` is a CeleScope pipeline for single-cell SNP/variant detection. It combines scRNA-Seq data with targeted gene panels to identify genetic variants at single-cell resolution.

This pipeline is designed to work with targeted capture technologies such as [FocuSCOPE](https://singleron.bio/products/focuscope-single-cell-targeted-capture-kit/), a single-cell targeted capture kit that enriches specific gene regions of interest while retaining full transcriptome information.

The pipeline runs the following steps in order:

1. **sample** — inspect FASTQ files and detect chemistry.
2. **starsolo** — align reads with STARsolo and generate expression matrices.
3. **target_metrics** — filter BAM to cell-associated reads mapped to target genes, and perform UMI deduplication.
4. **variant_calling** — call variants with bcftools after GATK SplitNCigarReads.
5. **filter_snp** — filter low-confidence alleles and update genotypes.
6. **analysis_snp** — annotate variants with snpEff and generate genotype matrices.

---

## Usage

### Make a snp reference genomeDir

1. Run `celescope rna mkref`. If you already have a rna genomeDir, you can use it and skip this step.
2. Run `celescope snp mkref` under the rna genomeDir.  This command will create [dictionary file and fasta index]((https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) ) for gatk SplitNCigarReads.
```
# under rna genomeDir
celescope snp mkref
```

### Run multi_snp

```
multi_snp\
    --mapfile ./test1.mapfile\
    --genomeDir {genomeDir after running celescope snp mkref}\
    --thread 16\
    --mod shell\
    --panel lung_1\
```

`--mapfile` Required.  Mapfile is a tab-delimited text file with for columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files. 
4th column: The single cell rna directory after running CeleScope.

Example

```tsv
R230630009      ./fastqs        TCR_sample_name  some_path/celescope_scrna_root/sample_name/
```

**Gene Selection Parameters & Priority**

--panel

    Description: Selects a predefined target gene panel (shortcut). The pipeline maps the panel name to its corresponding internal BED file to extract the gene list for downstream filtering and analysis.

    Available Options: e.g., lung_1, blood_1, CHIP.

    Use Case: Ideal when using standard gene sets without providing manual files.

--gene_list

    Description: A user-provided text file containing a list of gene symbols (one gene per line).

    Function: The program reads this file directly to define the gene set used for filtering, annotation, and variant retention.

--bed

    Description: Directly specifies a BED file containing genomic coordinates of genes or regions of interest.

    Function: The program extracts the gene names/regions directly from the BED file.

Logic and Impact
Priority Level

The pipeline selects the gene source based on the following order of precedence:

    --bed (Highest priority)

    --panel (Used if no BED is provided)

    --gene_list (Used if neither BED nor Panel is provided)

    Note: At least one of these three parameters must be provided, or the program will return an error.


**snpeff database**

We have identified an issue regarding variant annotation where SnpEff fails to download or access databases due to a change in their official repository URLs. This is a known issue discussed in the SnpEff GitHub community [Issue#602](https://github.com/pcingola/SnpEff/issues/602).

To resolve this and ensure the analysis_snp step runs correctly, please use one of the following two solutions:

Solution 1: Update the snpEff.config file (Recommended)

You can globally update the database URL in your environment's configuration file using a sed command. This ensures all future runs use the correct server.

Run the following command:

```Bash

sed -i 's|https://snpeff.blob.core.windows.net|https://snpeff.odsp.astrazeneca.com|g' /SGRNJ01/Public/Software/conda_env/celescope/share/snpeff-5.2-1/snpEff.config
```

Note: Please modify the path to your snpEff.config.

Solution 2: Use the `--snpeff_cache` parameter

If you prefer not to modify the global configuration or do not have write permissions, you can manually download the required database and specify the path during the multi_snp execution.

Action:
Add the `--snpeff_cache` argument to your multi_snp command to point to your local database directory:

```Bash
multi_snp \
    --snpeff_cache /path/to/your/local/snpeff_data \
    ... [other parameters]
```

The snpeff_cache for `GRCh38.mane.1.0.ensembl` is available for download at:
https://singleronbio-opendata.oss-cn-hangzhou.aliyuncs.com/genome/snpeff/snpeff_cache.tar


## Steps

### STARsolo (Mapping & Demultiplexing)
#### Purpose

Standard preprocessing of raw sequencing data to identify cell origins and gene targets.
#### Process

    Aligns raw FASTQ files to the reference genome.

    Assigns reads to individual cells using Barcodes and UMIs.

    Performs initial cell calling and filtering.

#### Output

    Gene expression matrices: Filtered matrices and raw count data.

    Primary mapping statistics.

### Target Metrics (Filtering & Enrichment)
#### Purpose

Refines the data by focusing on specific regions of interest and removing technical noise.
#### Process

    Filters the BAM file to retain only reads matching valid Cell Barcodes and target genes (via --gene_list or --panel).

    Performs UMI deduplication based on the combination of (CB, UMI, Reference_Position).

    Retains a maximum number of duplicates per UMI as specified.

#### Output

    {sample}_filtered.bam: A filtered and indexed BAM file for variant calling.

    {sample}_gene_count.tsv: Gene-level enrichment metrics.

### Variant Calling (SNP Identification)
#### Purpose

Identifies genetic variations (SNPs/Indels) directly from the processed RNA-seq data.
#### Process

    Applies GATK SplitNCigarReads to correctly handle RNA-seq spliced alignments.

    Calls variants using bcftools mpileup followed by bcftools call.

    Normalizes and left-aligns variants for consistent nomenclature.

#### Output

    {sample}_raw.vcf: The raw variant call file.

    {sample}_norm.vcf: A normalized VCF file with left-aligned variants.

### Filter SNP (Quality Control)
#### Purpose

Ensures high confidence in cell-specific genotypes by removing low-quality calls.

#### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bcftools_filter` | bcftools filter expression applied before allele-level filtering. | `QUAL>=100` |
| `--ref_threshold_method` | Threshold method for reference allele: `otsu` or `none`. | `none` |
| `--alt_threshold_method` | Threshold method for alternate allele: `otsu` or `none`. | `none` |
| `--ref_min_support_read` | Minimum supporting reads for reference allele. | `2` |
| `--alt_min_support_read` | Minimum supporting reads for alternate allele. | `2` |
| `--vaf` | VAF threshold for genotype calling. If `vaf <= ALT/(REF+ALT) <= 1-vaf`, genotype is heterozygous (`0/1`); otherwise homozygous (`0/0` or `1/1`). | `0.2` |

#### Algorithm

1. **Variant-level filtering**: Apply `bcftools filter -i` with `--bcftools_filter` to remove low-quality variant sites.

2. **Allele-level filtering**: For each variant site, collect reference and alternate allele depths (AD) across all cells. Determine a threshold for each allele using the specified method:
   - **`none`**: No statistical thresholding; only `--min_support_read` is applied.
   - **`otsu`**: Apply Otsu's method on the log-transformed positive allele counts to separate signal from noise. The threshold is computed per variant and per allele. When `otsu` is selected, histogram plots are saved to `ref_otsu_plots/` and `alt_otsu_plots/`.

   The final threshold is the maximum of the statistical threshold and `--min_support_read`.

3. **Count filtering**: Set AD counts below the final threshold to zero.

4. **Genotype calling**: Update genotypes (GT) based on the filtered AD counts:

   | Condition | Result (GT) |
   |-----------|-------------|
   | REF = 0, ALT = 0 | `NA` |
   | REF > 0, ALT = 0 | `0/0` |
   | REF = 0, ALT > 0 | `1/1` |
   | REF > 0, ALT > 0, `vaf <= AF <= 1 - vaf` | `0/1` |
   | REF > 0, ALT > 0, `AF < vaf` | `0/0` |
   | REF > 0, ALT > 0, `AF > 1 - vaf` | `1/1` |

   `AF = ALT / (REF + ALT)` is the alternate allele frequency after filtering.

#### Output

    {sample}_filtered.vcf: A filtered VCF file with updated AD and GT fields.

### Analysis SNP (Annotation & Visualization)
#### Purpose

Provides biological context and generates interpretable results for downstream analysis.
#### Process

    Annotates variants (e.g., missense, synonymous) using snpEff.

    Generates a Genotype Matrix mapping variants to individual cells.


#### Output

    {sample}_gt.csv: A Genotype Matrix mapping variants to individual cells.

    {sample}_variant_table.csv: Annotated variant table with snpEff information.