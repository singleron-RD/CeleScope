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

getwd()
print(matrix_dir)
mtx = Seurat::Read10X(matrix_dir)
bcrank = barcodeRanks(mtx)
uniq <- !duplicated(bcrank$rank)
res$FDR[is.na(res$FDR)] = 1
res.valid = res[res$FDR >0 & res$FDR<0.01,]
res.valid.ncell = dim(res.valid)[1]
print(paste('rescued cell number: ',res.valid.ncell))

res.out = stringr::str_glue('{outdir}/{sample}_rescue.tsv')
write.table(res, res.out, sep='\t')
res.valid.out = stringr::str_glue('{outdir}/{sample}_rescue_valid.tsv')
write.table(res.valid, res.valid.out, sep='\t')




