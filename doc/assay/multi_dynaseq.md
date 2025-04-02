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

- `{sample}_labeled_detail.csv`  comma-delimited file:
    - Barcode: Cell barcode sequence
    - UMI: UMI sequence
    - geneID: gene ID
    - TC: TC site number in a read (backgroup snp removed)


## Arguments

`--genomeDir` Required. Genome directory after running `celescope {assay} mkref`.

`--basequalilty` Min base quality of the read sequence. (20)

`--snp_threshold` SNP threshold filter, greater than snp_threshold will be recognized as snp. (0.5)

`--snp_min_depth` Minimum depth to call a variant. (20)

`--readsplit` Split to have approximately N reads per output file. (1,000,000)

`--all_type_plot` Plot subsititution rate for all conversion type.

`--control` For control samples to generate backgroup snp files and skip replacement.

`--conversionMem` Set conversion memory (G). (5)

`--replacementMem` Set replacement memory (G). (50)
