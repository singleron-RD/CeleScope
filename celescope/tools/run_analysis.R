library(Seurat) # v4.0
library(tidyverse)
library(argparser)
library(hdf5r)
library(rhdf5)


argv <- arg_parser('')
argv <- add_argument(argv,"--matrix_file", help="matrix file")
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- add_argument(argv,"--sample", help="sample")
argv <- add_argument(argv,"--save_rds", help="write rds to disk")
argv <- parse_args(argv)

#args
matrix_file = argv$matrix_file
outdir = argv$outdir
sample = argv$sample
save_rds = argv$save_rds
resolution = 0.6
res_str = paste0('res.', resolution)

tsne.out = stringr::str_glue('{outdir}/{sample}_tsne_coord.tsv')
marker.out = stringr::str_glue('{outdir}/{sample}_markers.tsv')
mito.out = paste(outdir,"stat.txt",sep="/")
rds.out = paste0(outdir,'/',sample,'.rds')

# read 10X
matrix = read.table(matrix_file,sep="\t",header=TRUE,row.names=1,quote = "")
rds = CreateSeuratObject(matrix, pro=sample)

# generate h5ad file
x = GetAssayData(rds,slot="count")
mtx = as.matrix(x)
barcode = colnames(rds)
geneid = rownames(rds)
h5.out = stringr::str_glue('{outdir}/{sample}.h5')
path <- path.expand(h5.out)
h5createFile(path)
h5f <- H5Fopen(path)
h5writeDataset(mtx,h5f,"X")
h5writeDataset(barcode,h5f,"obs")
h5writeDataset(geneid,h5f,"var")
H5Fclose(h5f)

# mito
mito.genes <- grep(pattern = "^MT-", x = rownames(x = rds@assays$RNA@data), value = TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(rds@assays$RNA@counts[mito.genes,])/Matrix::colSums(rds@assays$RNA@counts)
rds <- AddMetaData(object = rds, metadata = percent.mito, col.name = "percent.mito")
meta = rds@meta.data
total_cell = dim(meta)[1]
percent_list = c(0.05,0.1,0.15,0.2,0.5)
mito_df = dplyr::tibble(mito_percent=numeric(),cell_percent=numeric())
for (percent in percent_list){
  cell_percent = sum(meta$percent.mito > percent) / total_cell
  mito_df = mito_df %>% dplyr::add_row(mito_percent=percent,cell_percent=cell_percent)
}
paste0(round(mito_df$cell_percent * 100,2),"%")
mito_df$cell_percent = paste0(round(mito_df$cell_percent * 100,2),"%")
mito_df$mito_percent = paste0("Fraction of cells have mito gene percent>",round(mito_df$mito_percent * 100,2),"%")
write_delim(mito_df, mito.out, col_names=F, delim=":")


rds <- NormalizeData(rds, normalization.method = "LogNormalize",scale.factor = 10000)
rds <- FindVariableFeatures(rds, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf),
                            mean.function = ExpMean, dispersion.function = LogVMR)

use.genes <- rds@assays$RNA@var.features
rds <- ScaleData(rds, vars.to.regress = c("nCount_RNA", "percent.mito"), features = use.genes)
rds <- RunPCA(object = rds, features = use.genes, do.print = FALSE)
rds <- FindNeighbors(rds, dims = 1:20, force.recalc = TRUE, reduction = "pca")
rds <- FindClusters(rds, resolution = resolution)

# tsne and umap
rds <- RunTSNE(rds, dims = 1:20, do.fast = TRUE, check_duplicates = FALSE)


tryCatch({
  rds.markers <- FindAllMarkers(object = rds, features = use.genes)
  rds.markers = dplyr::group_by(rds.markers,cluster) %>% dplyr::arrange(desc(avg_log2FC))
}, error = function(e){
  print (paste0("no marker found: ", e))
  rds.markers <<- data.frame(cluster=double(),
                             gene=double(),
                             avg_log2FC=double(),
                             pct.1=double(),
                             pct.2=double(),
                             p_val_adj=double())

})

rds.markers$cluster = as.numeric(rds.markers$cluster)
print (rds.markers)
write_tsv(rds.markers,marker.out,col_names = T)


df.tsne = rds@reductions$tsne@cell.embeddings
df.tsne = as.data.frame(df.tsne)
meta = rds@meta.data
dic = rds@meta.data[['seurat_clusters']]
names(dic) = rownames(rds@meta.data)
df.tsne$cluster = as.numeric(dic[rownames(df.tsne)])
rds@meta.data$seurat_clusters = as.numeric(dic[rownames(df.tsne)])
df.gene = meta[,"nFeature_RNA",drop=F]
colnames(df.gene) = "Gene_Counts"
df.all = cbind(df.tsne,df.gene)
write.table(df.all,tsne.out,sep="\t",col.names=NA,quote = F)


if (save_rds == 'True'){
  saveRDS(rds, rds.out)
}