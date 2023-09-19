## Usage

```
    multi_dynaseq\\
    --mapfile ./rna.mapfile\\
    --genomeDir /SGRNJ/Public/Database/genome/homo_mus
```

For control sample, set --control to skip replacement step.
```
    multi_dynaseq\\
    --mapfile ./rna.mapfile\\
    --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
    --control
```

For genome reference generated, please refer to rna assay (./multi_rna.md) .

## Main Output
- `substitution/{sample}.substitution.txt`: Substitution rate for each conversion type.
- `replacement/{sample}.labeled.h5ad`: h5ad file contains ['total', 'labeled', 'unlabeled'] layers.

## Features
### featureCounts
- Assigning uniquely mapped reads to genomic features with FeatureCounts.

### count
- Cell-calling: Distinguish cell barcodes from background barcodes. 
- Generate expression matrix.

### analysis
- Cell clustering with Seurat.
- Calculate the marker gene of each cluster.
- Cell type annotation(optional). You can provide markers of known cell types and annotate cell types for each cluster.

### conversion
- Get conversion pos in each read.
    - Get snp info. 

### substitution
- Computes the overall conversion rates in reads and plots a barplot.

### replacement
- Quantify unlabeled and labeled RNA.
- Boxplots for TOR rates distribution.
- TSNE plot for TOR rate 


## Output files
### featureCounts
- `featureCounts/{sample}_aligned_posSorted_addTag.bam` This bam file contains coordinate-sorted reads aligned to the genome.

### count
- `{sample}_raw_feature_bc_matrix` The expression matrix of all detected barcodes in [Matrix Market Exchange Formats](
    https://math.nist.gov/MatrixMarket/formats.html). 
- `{sample}_filtered_feature_bc_matrix` The expression matrix of cell barcodes in Matrix Market Exchange Formats. 

### analysis
- `markers.tsv` Marker genes of each cluster.

### conversion
- `{sample}.PosTag.bam` Bam file with conversion info.
- `{sample}.PosTag.csv` TC conversion sites info in csv format.
- `{sample}.snp.csv` Candidated snp sites.

### substitution
- `{sample}.substitution.txt` Tab-separated table of the overall conversion rates.

### replacement
- `{sample}.labeled.h5ad` h5ad file contains ['total', 'labeled', 'unlabeled'] layers and TOR rate of each cell/gene.
- `{sample}_labeled_feature_bc_matrix` The labeled expression matrix of cell barcodes & all features in Matrix Market Exchange Formats. (will be deprecated in future versions)
- `{sample}_unlabeled_feature_bc_matrix` The unlabeled expression matrix of cell barcodes & all features in Matrix Market Exchange Formats. (will be deprecated in future versions)
- `{sample}_labeled_detail.txt`  tab-delimited  file:
    - Barcode: Cell barcode sequence
    - UMI: UMI sequence
    - geneID: gene ID
    - TC: TC site number in a read (backgroup snp removed)

## Arguments

`--genomeDir` Required. Genome directory after running `celescope {assay} mkref`.

`--basequalilty` Min base quality of the read sequence.

`--snp_min_cells` Minimum number of cells to call a variant(>=1 for cell number or <1 for cell fraction).

`--snp_min_depth` Minimum depth to call a variant.

`--cellsplit` Split N cells into a list.

`--conversionMem` Set conversion memory.

`--control` For control samples to generate backgroup snp files and skip replacement.

