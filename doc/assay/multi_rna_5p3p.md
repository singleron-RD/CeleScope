Before using, please refer to [multi_rna.md](./multi_rna.md) to learn how to create genomeDir and mapfile.

Using `multi_rna_5p3p` is very similar to `multi_rna`, with the following differences:
- Each sample has at least two FASTQ files: 5-prime and 3-prime, so the mapfile requires a fourth column(`5p` or `3p`) to specify whether the FASTQ file is 5-prime or 3-prime.
- You must manually specify `--chemistry 5p3p-1`.

**mapfile Format**
```tsv
fastq_prefix_5p fastq_5p_folder sample_name 5p
fastq_prefix_3p fastq_3p_folder sample_name 3p
```

**run.sh Example**
```shell
multi_rna_5p3p \
 --mapfile 5p3p.mapfile \
 --genomeDir path_of_genomeDir \
 --chemistry 5p3p-1 \
 --mod shell
```


## Downstream analysis using [MARVEL](https://wenweixiong.github.io/MARVEL_Droplet.html)
After running celescope rna_5p3p, the output directory will contain the following directory structure:
```
00.sample
01.convert
02.starsolo
03.analysis
outs
```

### Required Files
#### Gene Matrix
This can be found at `outs/filtered`.

#### SJ Matrix
The SJ matrix is located at `02.starsolo/*_Solo.out/SJ/raw`. In order for the SJ matrix's features file to be compatible with MARVEL analysis, you need to process the raw features.tsv file. You can use the following command:
```
awk '{print "chr"$1":"$2":"$3}' SJ_raw_features.tsv > SJ_raw_features.tsv
```
This command modifies the feature names in the original SJ matrix features.tsv to the format chromosome:intron_start_position:intron_end_position, with the chromosome name prefixing as chr. The chromosome names in the gtf file must also have the chr prefix.

#### GTF File
Use the corresponding genome's GTF file from the celescope rna_5p3p pipeline. The chromosome names in the GTF file should have the chr prefix and should not include special characters like `-`.
#### SJ Read Count File
This can be found at `02.starsolo/*_SJ.out.tab`.

### Optional Files
#### Dimensionality Reduction Coordinates File
This file is located at `outs/tsne_coord.tsv`. If this file is not provided, MARVEL can use Seurat to re-cluster the data and generate new coordinate files.