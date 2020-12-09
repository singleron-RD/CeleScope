library(DropletUtils)
library(Seurat)
library(argparser)
library(tidyverse)

argv <- arg_parser('')
argv <- add_argument(argv,"--outdir", help="the output dir.",default=getwd())
argv <- add_argument(argv,"--sample", help="sample")
argv <- add_argument(argv,"--matrix_dir", help="matrix_dir")
argv <- add_argument(argv,"--threshold", help="UMI threshold")
argv <- parse_args(argv)

#read args
outdir = argv$outdir
sample = argv$sample
matrix_dir = argv$matrix_dir
threshold = as.numeric(argv$threshold)

mtx = Seurat::Read10X(matrix_dir)
res = emptyDrops(mtx, lower=threshold, niters=1000, test.ambient=TRUE, retain=threshold)
res$FDR[is.na(res$FDR)] = 1
res = res[res$FDR >0 & res$FDR<0.01,]
res.ncell = dim(res)[1]
print(paste('rescued cell number: ',res.ncell))
res.out = stringr::str_glue('{outdir}/{sample}_rescue.tsv')
write.table(res, res.out)




