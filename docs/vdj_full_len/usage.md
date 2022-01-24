# Example

Test path: 

`/SGRNJ03/randd/cjj/celedev/TESTDATA/test_cellranger`

Conda environment: 

`cellranger_VDJ`

Run command: 

```
multi_vdj_full_len \
    --mapfile  ./test.mapfile \
    --chemistry flv \
    --mem 10\
    --species hs \
    --thread 8 \
    --allowNoLinker \
    --soft 4.0.0 --seqtype TCR --mod sjm
```

```
--chemistry flv
--species Choose from ['hs', 'mmu']. Required.
--soft Select soft version. chose from ['3.0.2', '3.1.0', '4.0.0', '6.0.0']. Default 4.0.0.
--mem Memory(G) for assembly. Default 10.
--thread number of multi-threaded for assembly.
--seqtype Choose from ['TCR', 'BCR']. Required.
```
