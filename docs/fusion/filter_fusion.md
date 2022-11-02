## Features
- Correct single-base errors in UMIs due to sequencing, amplification, etc.
- Filter background UMIs base on a UMI threshold.
There are three methods to determine the UMI threshold:
    - 'auto' : Using a method similar to cell calling method.
    - 'otsu' : UMI counts are first log 2 transformed and then the threshold is determined by [Otsu's method](https://en.wikipedia.org/wiki/Otsu%27s_method)
    - 'hard' : Using User provided UMI threshold.

## Output
- `{sample}_corrected_read_count.json` Read counts after UMI correction.
- `{sample}_filtered_read_count.json` Filtered read counts.
- `{sample}_filtered_UMI.csv` Filtered UMI counts.
## Arguments
`--not_correct_UMI` Do not perform UMI correction.

`--read_threshold_method` method to find read threshold. UMIs with `support reads` < `read threshold` are filtered.

`--read_hard_threshold` int, use together with `--read_threshold_method hard`.

`--umi_threshold_method` method to find UMI threshold. Cell barcode with `UMI` < `UMI threshold` are considered negative.

`--umi_hard_threshold` int, use together with `--umi_threshold_method hard`.

`--auto_coef` int, threshold = top 1 percent positive cell count / auto_coef.

`--otsu_log_base` raw counts are first log transformed before thresholding. This argument is the log base. Commonly used values are 2 and 10.

`--match_dir` Match celescope scRNA-Seq directory.

`--raw_read_count_file` None

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

