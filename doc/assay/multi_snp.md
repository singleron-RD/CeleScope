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



