# Plot variants genotype on t-SNE

## Usage

```
Rscript plot_snp.R \
 -d /SGRNJ03/randd/user/zhouyiqi/tests/19del_snp/Cap211019a_GC_FJ \
 -o ./output_dir \
 --protein 692_697del, T790M \
 --position 7_55174771, 7_55181378
```

`-d, --snp_dir` CeleScope snp directory.

`-o, --outdir` Output directory.

`--protein` Optional. Plot variants with these protein changes. Multiple protein changes are separated by comma.

`--position` Optional. Plot variants at these positions. Multiple positions are separated by comma.

## Output

t-SNE plot with genotypes(NA, 0/0, 0/1, 1/1). 