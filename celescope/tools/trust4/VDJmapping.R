library(Seurat)
library(tidyverse)
library(argparse)
library(ggplot2)
library(stringr)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--sample",help="sample name", required=TRUE)
parser$add_argument("--rds", help="scRNA rds path", required=TRUE)
parser$add_argument("--VDJ", help='VDJ match_contigs.csv path', required=TRUE)
parser$add_argument("--outdir", help='out dir', required=TRUE)
parser$add_argument("--assign_file", help='auto_assign', required=TRUE)
args <- parser$parse_args()

rna <- readRDS(args$rds)
vdj <- read.table(args$VDJ, sep=',', header=T)

cells <- subset(vdj, productive=='True')
barcodes <- unique(cells$barcode)

df <- rna@meta.data
df$barcode <- rownames(df)

filter_df <- filter(df, barcode %in% barcodes)
res <- table(filter_df$seurat_clusters)
res <- as.data.frame(res)

out_csv = stringr::str_glue("{args$outdir}/{args$sample}_match_barcodes_celltypes_distribution.txt")
write_csv(res,file=out_csv)

meta = rna@meta.data
meta$Class = 'NA'
meta[barcodes,'Class'] = 'T/BCR'
rna@meta.data = meta
rna <- RunUMAP(rna, dims = 1:20)

outP = stringr::str_glue("{args$outdir}/{args$sample}_cluster_umap.png")
png(outP, height=1000, width=1000)
UMAPPlot(rna,group.by='seurat_clusters',label=TRUE)
dev.off()

outP1 = stringr::str_glue("{args$outdir}/{args$sample}_umapplot.png")
png(outP1, height=1000, width=1000)
UMAPPlot(rna,group.by='Class',cols=c('grey','red'),label=TRUE)
dev.off()

cell_ident_file <- read.table(args$assign_file,header = TRUE,sep="\t",stringsAsFactors=FALSE)
cell_ident_file <- cell_ident_file %>%
  group_by(cluster) %>%
  filter(avg_pct.diff==max(avg_pct.diff))
current_ident <- cell_ident_file[,1]
new_ident <- cell_ident_file[,2]
new.cluster.ids = (new_ident$cell_type)
Idents(rna) <- "seurat_clusters"
names(new.cluster.ids) <- sort(as.numeric(levels(rna)))
rna <- RenameIdents(rna, new.cluster.ids)
rna <- StashIdent(object = rna, save.name = "CellTypes")

outP2 = stringr::str_glue("{args$outdir}/{args$sample}_assign.png")
png(outP2, height=1000, width=1000)
UMAPPlot(rna,group.by='CellTypes',label=TRUE,label.box=TRUE)
dev.off()

meta = rna@meta.data
outMeta = stringr::str_glue("{args$outdir}/{args$sample}_meta.csv")
write.csv(meta,outMeta)

