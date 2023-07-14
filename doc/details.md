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
Please note that only cell barcodes are considered and all the background barcodes are discarded. So the cell-calling results on each tag may be different from the results obtained from `--split_matrix` and `--split_vdj`.

1.  Run `barcode` to get demultiplexed R2 fastq. If you have already run `multi_rna`, `multi_vdj`, `multi_dynaseq`, etc.., this step can be skipped.
```
multi_rna \ 
 --mapfile ./rna.mapfile \
 --mod shell \
 --steps_run sample,barcode \
 --genomeDir {genomeDir}
```

2. Run `multi_tag`
```
multi_tag \
 --mapfile ./tag.mapfile\
 --mod shell\
 --barcode_fasta ./tag_barcode.fasta\
 --fq_pattern L25C15 \
 --split_fastq \
```

3. Add `--R1_read` in `./shell/{sample}.sh`
```
...
...
celescope tag analysis_tag ...
celescope tag split_tag ... --R1_read {e.g. R1_read of the matched scRNA-Seq library}
```
The output fastqs are in `{sample}/06.split_tag/`


