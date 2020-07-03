library(Seurat)
library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--matrix_file", help="matrix file")
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- add_argument(argv,"--sample", help="sample")
argv <- parse_args(argv)

#args
matrix_file = argv$matrix_file
outdir = argv$outdir
sample = argv$sample

matrix = read.table(matrix_file,sep="\t",header=TRUE,row.names=1)
tsne.out = paste(outdir,"tsne_coord.tsv",sep="/")
marker.out = paste(outdir,"markers.tsv",sep="/")
mito.out = paste(outdir,"stat.txt",sep="/")


rds = CreateSeuratObject(raw.data = matrix,project=sample)

# mito
mito.genes <- grep(pattern = "^MT-", x = rownames(x = rds@data), value = TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(rds@raw.data[mito.genes,])/Matrix::colSums(rds@raw.data)
rds <- AddMetaData(object = rds, metadata = percent.mito, col.name = "percent.mito")
meta = rds@meta.data
total_cell = dim(meta)[1]
percent_list = c(0.05,0.1,0.15,0.2)
mito_df = dplyr::tibble(mito_percent=numeric(),cell_percent=numeric())
for (percent in percent_list){
  cell_percent = sum(meta$percent.mito > percent) / total_cell
  mito_df = mito_df %>% dplyr::add_row(mito_percent=percent,cell_percent=cell_percent)
}
paste0(round(mito_df$cell_percent * 100,2),"%")
mito_df$cell_percent = paste0(round(mito_df$cell_percent * 100,2),"%")
mito_df$mito_percent = paste0("Fraction of cells have mito gene percent>",round(mito_df$mito_percent * 100,2),"%")
write_tsv(mito_df,"mito.tsv",col_names = F)

rds <- NormalizeData(object = rds, normalization.method = "LogNormalize",scale.factor = 10000)
rds <- FindVariableGenes(object = rds, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1,  y.cutoff = 1)
use.gene = rds@var.genes
rds <- ScaleData(object = rds,vars.to.regress = c("nUMI", "percent.mito"),genes.use =use.gene)
rds <- RunPCA(object = rds, pc.genes = use.gene, do.print = FALSE)
rds <- FindClusters(object = rds, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE,force.recalc = TRUE)

# Run Non-linear dimensional reduction (tSNE)
rds <- RunTSNE(object = rds, dims.use = 1:20, do.fast = TRUE,check_duplicates = FALSE)
rds.markers <- FindAllMarkers(object = rds, genes.use = use.gene)
rds.markers = dplyr::group_by(rds.markers,cluster) %>% dplyr::arrange(desc(avg_logFC))

write_tsv(rds.markers,marker.out)

df.tsne = rds@dr$tsne@cell.embeddings
df.tsne = as.data.frame(df.tsne)
meta = rds@meta.data
dic = rds@meta.data$res.0.6
names(dic) = rownames(rds@meta.data)
df.tsne$cluster = as.character(dic[rownames(df.tsne)])
df.gene = meta[,"nGene",drop=F]
colnames(df.gene) = "Gene_Counts"
df.all = cbind(df.tsne,df.gene)
write.table(df.all,tsne.out,sep="\t",col.names=NA,quote = F)