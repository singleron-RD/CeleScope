# General
## Valid reads

As you can see from https://github.com/singleron-RD/CeleScope/blob/master/docs/tools/barcode.md,  valid reads are based on R1 reads. The reads are filtered in the following order: 

`no polyT - low quality - incorrect linker - incorrect barcode.` 

If a read passes all four filtering steps, it is classified as a valid read.

You can check the log of the `barcode` step for the number of filtered reads in each filtering step.
```
2022-03-29 09:58:47,733 - celescope.tools.barcode.run - INFO - processed reads: 100,000. valid reads: 73,300.
2022-03-29 09:58:47,733 - celescope.tools.barcode.run - INFO - no polyT reads number : 23196
2022-03-29 09:58:47,733 - celescope.tools.barcode.run - INFO - low qual reads number: 0
2022-03-29 09:58:47,733 - celescope.tools.barcode.run - INFO - no_linker: 2431
2022-03-29 09:58:47,734 - celescope.tools.barcode.run - INFO - no_barcode: 1073
2022-03-29 09:58:47,734 - celescope.tools.barcode.run - INFO - corrected linker: 5642
2022-03-29 09:58:47,734 - celescope.tools.barcode.run - INFO - corrected barcode: 4728
```

If `no polyT reads number` is very high, use `--allowNoPolyT`. This may be caused by excessive sequencing errors in the polyT region.

PS: If the default parameters are used, low quality reads are not filtered. So `low qual reads number` will always be zero.

## Chemistry

### The correspondence between chemistry and kits

- `scopeV1` : Micro Bead Kit. This kit is no longer in use.
- `scopeV2.*.*` : Magnetic Bead Kit V1.
- `scopeV3.*.*` : Magnetic Bead Kit V2.

### Analyze data with the updated chemistry version: scopeV3.0.1
- If you are using the latest CeleScope version(>= v1.5.2), it can auto-detect scopeV3.0.1, so `--chemistry auto` will work. You can also explicitly specify `--chemistry scopeV3.0.1`.

- If you are using CeleScope version < v1.5.2, then you need to use `--chemistry customized` and provide 3 additional arguments: (pattern, whitelist, linker).
```
multi_rna \
 --chemistry customized \
 --pattern C9L16C9L16C9L1U12T18 \
 --whitelist celescope/data/chemistry/scopeV3.0.1/bclist \
 --linker celescope/data/chemistry/scopeV3.0.1/linker_4types \
 ...
```

- There are no changes in other arguments.

### Where to find pattern, whitelist and linker of each chemistry?
- pattern: https://github.com/singleron-RD/CeleScope/blob/master/celescope/tools/__init__.py
- whitelist and linker:  https://github.com/singleron-RD/CeleScope/tree/master/celescope/data/chemistry/ . `bclist` is the barcode whitelist.

## Saturation

There is a difference in how CeleScope and CellRanger calculate saturation. CeleScope shows umi_saturation in the report, while CellRanger shows read_saturation in the report. 

```
umi_saturation = 1 - (n_deduped_reads / n_umis)
read_saturation = 1 - (n_deduped_reads / n_reads)

n_deduped_reads = Number of unique (valid cell-barcode, valid UMI, gene) combinations among confidently mapped reads.

n_umis = Total number of (confidently mapped, valid cell-barcode, valid UMI) UMIs.

n_reads = Total number of (confidently mapped, valid cell-barcode, valid UMI) reads.
```

You can find the saturation calculated in both ways in {sample}/.metrics.json

```
{
    "count_summary": {
        "Read Fraction 0.1 Read_saturation": 2.62,
        "Read Fraction 0.1 Umi_saturation": 1.33,
        "Read Fraction 0.2 Read_saturation": 6.81,
        "Read Fraction 0.2 Umi_saturation": 3.47,
        "Read Fraction 0.3 Read_saturation": 10.61,
        "Read Fraction 0.3 Umi_saturation": 5.47,
        "Read Fraction 0.4 Read_saturation": 13.47,
        "Read Fraction 0.4 Umi_saturation": 7.0,
        "Read Fraction 0.5 Read_saturation": 16.05,
        "Read Fraction 0.5 Umi_saturation": 8.38,
        "Read Fraction 0.6 Read_saturation": 18.99,
        "Read Fraction 0.6 Umi_saturation": 10.06,
        "Read Fraction 0.7 Read_saturation": 21.52,
        "Read Fraction 0.7 Umi_saturation": 11.46,
        "Read Fraction 0.8 Read_saturation": 24.13,
        "Read Fraction 0.8 Umi_saturation": 12.99,
        "Read Fraction 0.9 Read_saturation": 26.27,
        "Read Fraction 0.9 Umi_saturation": 14.2,
        "Read Fraction 1.0 Read_saturation": 29.27,
        "Read Fraction 1.0 Umi_saturation": 16.07,
        "Estimated Number of Cells": 167,
        "Fraction Reads in Cells": 63.9,
        "Median UMI per Cell": 37,
        "Total Genes": 3321,
        "Median Genes per Cell": 34,
        "Saturation": 16.07
    }
}
```




# Clindex

## Tag assignment

[tag assignment algorithm](https://github.com/singleron-RD/CeleScope/blob/master/methods/tag_algorithm.txt)

You can control the value of UMI_min and SNR_min by changing the arguments in `celescope tag count_tag`
(https://github.com/singleron-RD/CeleScope/blob/master/docs/tag/count_tag.md)

The value of `UMI_min` and `SNR_min` is currently displayed in the log file of `celescope tag count_tag` and will be added to the report in subsequent releases.

## Split matrix
If you want to split the expression matrix of match scRNA-Seq library，you need to add `--split_matrix`  in the `multi_tag`, e.g.
```
multi_tag \
 --mapfile ./tag.mapfile\
 --mod shell\
 --barcode_fasta ./tag_barcode.fasta\
 --fq_pattern L25C15 \
 --split_matrix
```
The output matrices are in `{sample}/06.split_tag`

## Split fastq
To split the R1 and R2 fastq files:

1.  Run `barcode` to get demultiplexed R2 fastq
```
multi_rna \ 
 --mapfile ./rna.mapfile \
 --mod shell \
 --steps_run sample,barcode \
 --genomeDir {genomeDir}
```

As the `sample` and `barcode` are the same among all assays, it is ok to use `multi_dynaseq` and other `multi_{assay}` to demultiplex the R2 fastq.

2. Run `multi_tag`
```
multi_tag \
 --mapfile ./tag.mapfile\
 --mod shell\
 --barcode_fasta ./tag_barcode.fasta\
 --fq_pattern L25C15 \
 --split_fastq \
 --R1_read {R1 fastq path}
```
The output fastqs are in `{sample}/06.split_tag/`

