## Usage

```
multi_dynaseq \
    --mapfile ./rna.mapfile \
    --genomeDir /SGRNJ/Public/Database/genome/homo_mus
```

For control sample, set --control to skip replacement step.
```
multi_dynaseq \
    --mapfile ./rna.mapfile \
    --genomeDir /SGRNJ/Public/Database/genome/homo_mus \
    --control
```

For genome reference generated, please refer to [rna](multi_rna.md) assay.

## Main Output
- `outs/labeled` labeled matrix
- `outs/unlabeled` unlabeled matrix
- `outs/{sample}.labeled.h5ad`: h5ad file contains ['total', 'labeled', 'unlabeled'] layers.
- `substitution/{sample}.substitution.txt`: Substitution rate for each conversion type.

## Steps
### conversion
- Get conversion pos in each read.
    - Get snp info. 

### substitution
- Computes the overall conversion rates in reads and plots a barplot.

### replacement
- Quantify unlabeled and labeled RNA.
- Boxplots for TOR rates distribution.
- TSNE plot for TOR rate 


### conversion
- `{sample}.PosTag.bam` Bam file with conversion info.
- `{sample}.PosTag.csv` TC conversion sites info in csv format.
- `{sample}.snp.csv` Candidated snp sites.

### substitution
- `{sample}.substitution.txt` Tab-separated table of the overall conversion rates.

### replacement
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

