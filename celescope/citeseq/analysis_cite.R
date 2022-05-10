library(Seurat)
library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--citeseq_mtx", help="citeseq_mtx file")
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- add_argument(argv,"--sample", help="sample")
argv <- add_argument(argv,"--rds", help="rds")
argv <- parse_args(argv)

#args
citeseq_mtx = argv$citeseq_mtx
outdir = argv$outdir
sample = argv$sample
rds = argv$rds

rds = readRDS(rds)
rds.adt <- read.csv(citeseq_mtx, sep = "\t", header = TRUE, row.names = 1)
adt_assay <- CreateAssayObject(counts = rds.adt)
rds[["ADT"]] <- adt_assay

# Normalize ADT data,
DefaultAssay(rds) <- "ADT"
rds <- NormalizeData(rds, normalization.method = "CLR", margin = 2)
DefaultAssay(rds) <- "RNA"

# Note that the following command is an alternative but returns the same result
rds <- NormalizeData(rds, normalization.method = "CLR", margin = 2, assay = "ADT")

DefaultAssay(rds) <- "ADT"

pdf.feature = str_glue('{outdir}/{sample}_CITESeq_featurePlot.pdf')
tags = rownames(rds@assays$ADT@data)
n_tags = length(tags)

pdf(pdf.feature, height=n_tags, width=12)
p = FeaturePlot(rds, rownames(rds@assays$ADT@data), ncol=4, pt.size=0.2, reduction='tsne')
print(p)
dev.off()

pdf.vln = str_glue('{outdir}/{sample}_CITESeq_vlnPlot.pdf')
pdf(pdf.vln, height=n_tags, width=18)
p = VlnPlot(rds, rownames(rds@assays$ADT@data), ncol=4, pt.size=0.2)
print(p)
dev.off()