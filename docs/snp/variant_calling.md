## Features
- Perform variant calling at single cell level.

## Output

`{sample}_VID.tsv` A unique numeric ID is assigned for each variant.

`{sample}_CID.tsv` A unique numeric ID is assigned for each cell.

`{sample}_variant_ncell.tsv` VID count of ref and alt. `VID`: A unique numeric ID is assigned for each variant. `ncell_cover`: means number of cells with read count at this position. `ncell_alt`: means number of cells with variant read count only. `ncell_ref`: number of cells with reference read count only. `ncell_ref_and_alt`: means number of cells with meanwhile have read count and reference read count.

`{sample}_merged.vcf ` VCF file containing all variants of all cells. `VID` and `CID` are added to the `INFO` column.

`{sample}_filter.vcf` VCF file after filtering. Invalid `CID`s are removed from the `INFO` column.

`{sample}_variant_count.tsv`  Reference and variant supporting reads/UMIs counts.

`{sample}_filter_variant_count.tsv`  Reference and variant supporting reads/UMIs counts after filtering.

`{sample}_support.mtx` Support matrix in [Matrix Market Exchange Formats](https://math.nist.gov/MatrixMarket/formats.html). Rows 
are variants(VID) and columns are cells(CID). The value can be 1, 2 or 3.

1 : all reads/UMIs at the position support the ref allele.  
2 : all reads/UMIs at the position support the alt allele.  
3 : one or more reads/UMIs support both the alt and the ref allele.  


## Arguments
`--genomeDir` Required. Genome directory after running `mkref`.

`--min_support_read` Minimum number of reads support a variant. If `auto`(default), otsu method will be used to determine this value.

`--bam` Input BAM file from step `target_metrics`.

`--match_dir` Match celescope scRNA-Seq directory.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

