library(Seurat)
library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--outdir", help="the output dir.",default=getwd())
argv <- add_argument(argv,"--sample", help="sample")
argv <- add_argument(argv,"--umi_tag",help="umi_tag")
argv <- add_argument(argv,"--matrix_10X",help="matrix_10X")
argv <- parse_args(argv)

outdir = argv$outdir
sample = argv$sample
umi_tag = argv$umi_tag
matrix_10X = argv$matrix_10X

df = read.csv(umi_tag, sep='\t', row.names = 1)
counts = Seurat::Read10X(matrix_10X)

df = df[1:length(df)-1]
df = t(df)

rds = CreateSeuratObject(counts = counts)
rds[["HTO"]] <- CreateAssayObject(counts = df)

rds <- NormalizeData(rds, assay = "HTO", normalization.method = "CLR")
rds <- HTODemux(rds, assay = "HTO")

# output
global_class = data.frame(table(rds$HTO_classification.global))
global_class.out = str_glue("{outdir}/{sample}_seurat_global_class.tsv")
write_tsv(global_class, global_class.out, col_names = F)

tag_count = data.frame(table(rds$HTO_classification))
tag_count.out = str_glue("{outdir}/{sample}_seurat_tag_count.tsv")
write_tsv(tag_count, tag_count.out, col_names = F)

tag = data.frame(rds$HTO_classification)
tag.out = str_glue("{outdir}/{sample}_seurat_tag.tsv")
write_tsv(tag, tag.out, col_names = F)


