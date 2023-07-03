## The correspondence between chemistry and kits

- `scopeV1` : Micro Bead Kit. This kit is no longer in use.
- `scopeV2.*.*` : Magnetic Bead Kit V1.
- `scopeV3.*.*` : Magnetic Bead Kit V2.

## Analyze data with the updated chemistry version: scopeV3.0.1
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

## Where to find pattern, whitelist and linker of each chemistry?
- pattern: https://github.com/singleron-RD/CeleScope/blob/master/celescope/tools/__init__.py
- whitelist and linker:  https://github.com/singleron-RD/CeleScope/tree/master/celescope/data/chemistry/ . `bclist` is the barcode whitelist.
