## Arguments
`--genomeDir` Required. Genome directory after running `celescope {assay} mkref`.

`--outFilterMatchNmin` Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.

`--out_unmapped` Output unmapped reads.

`--STAR_param` Additional parameters for the called software. Need to be enclosed in quotation marks. For example, `--{software}_param "--param1 value1 --param2 value2"`.

`--outFilterMultimapNmax` Default `1`. How many places are allowed to match a read at most.

`--starMem` Default `30`. Maximum memory that STAR can use.

`--fq` Required. R2 fastq file.

`--consensus_fq` A indicator that the input fastq has been consensused.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

`--fusion_genomeDir` fusion gene STAR index genome directory.

