## Analyze data with the updated chemistry version: scopeV3.0.1
- If you are using the latest CeleScope version(>= v1.5.2), it can auto-detect scopeV3.0.1, so `--chemistry auto` wiil work. You can also explicitly specify `--chemistry scopeV3.0.1`.

- If you are using CeleScope version < v1.5.2, then you need to use `--chemistry customized` and provide 3 additional arguments: (pattern, whitelist, linker).

- There are no changes in other arguments.

## Where to find pattern, whitelist and linker of each chemistry?
- pattern: https://github.com/singleron-RD/CeleScope/blob/master/celescope/tools/__init__.py
- whitelist and linker:  https://github.com/singleron-RD/CeleScope/tree/master/celescope/data/chemistry/ . `bclist` is the barcode whitelist.