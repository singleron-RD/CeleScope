# Tests

1. Get test data and scripts.
```
mkdir root_dir
cd root_dir
git clone https://github.com/singleron-RD/celescope_test_data.git
git clone https://github.com/singleron-RD/celescope_test_script.git
```

2. Modify genome arguments
- `celescope_test_script/rna/run_shell.sh` Modify `--genomeDir`

3. Run `pytest`
```
Install pytest
>>> pip install pytest
Run all
>>> python -m pytest -s ./tests/test_multi.py --test_dir {root_dir/celescope_test_script}
Run some tests
>>> python -m pytest -s ./tests/test_multi.py --test_dir {root_dir/celescope_test_script} --assays tag,vdj
```