library(Seurat)
library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--matrix_dir", help="")
argv <- add_argument(argv,"--outdir", help="")
argv <- add_argument(argv,"--sample", help="")
argv <- add_argument(argv,"--mt_gene_list_file", help="")
argv <- parse_args(argv)

matrix_dir = argv$matrix_dir
outdir = argv$outdir
sample = argv$sample
mt_gene_list_file = argv$mt_gene_list_file

# out
df.out = str_glue("{outdir}/{sample}_mt_UMI.tsv")

mtx = Read10X(matrix_dir)
mt_gene_list = read.table(mt_gene_list_file)[,1]

gene_valid = rownames(mtx)
gene_intersect = intersect(gene_valid, mt_gene_list)
cells = dim(mtx)[2]
mean_UMI = sort(round(rowSums(mtx[gene_intersect,]) / cells,3), decreasing = T)
df = as.data.frame(mean_UMI)
write.table(df, df.out, sep='\t', col.names=NA, quote = F)
