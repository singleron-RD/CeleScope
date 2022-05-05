## Features
- Count the number of reads and umis that 
    1. originate from cell barcodes;
    2. align to the fusion site and include flanking sequences of a certain length(default 20bp) on both sides of the fusion site.
## Arguments
`--fusion_genomeDir` Fusion genome directory.

`--flanking_base` Number of bases flanking the fusion position.

`--min_query_length` Minimum query length.

`--match_dir` Match celescope scRNA-Seq directory.

`--capture_bam` None

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

