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


## Features

### conversion

Get conversion info in each read.

- `{sample}.PosTag.bam` Bam file with conversion info.
- `{sample}.PosTag.csv` TC conversion sites info in csv format.
- `{sample}.snp.csv` Candidated snp sites.


### substitution

Calculate the overall substitution rates in reads.

- `{sample}.TC_substitution.txt`: Substitution rate of labeled and backgroud sites.
- `{sample}.substitution.txt (optional)` Tab-separated table of the overall conversion rates.


### replacement

Quantify unlabeled and labeled RNA.

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

`--conversionMem` Set conversion memory.

`--all_type_plot` Plot subsititution rate for all conversion type.

`--cellsplit` Split N cells into a list for program parallelism.

`--control` For control samples to generate backgroup snp files and skip replacement.

