library(Seurat)
library(tidyverse)
library(argparser)

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

matrix = read.table(matrix_file,sep="\t",header=TRUE,row.names=1)
tsne.out = stringr::str_glue('{outdir}/{sample}_tsne_coord.tsv')
marker.out = stringr::str_glue('{outdir}/{sample}_markers.tsv')
mito.out = paste(outdir,"stat.txt",sep="/")
rds.out = paste0(outdir,'/',sample,'.rds')


rds = CreateSeuratObject(raw.data = matrix,project=sample)

# mito
mito.genes <- grep(pattern = "^MT-", x = rownames(x = rds@data), value = TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(rds@raw.data[mito.genes,])/Matrix::colSums(rds@raw.data)
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

rds <- NormalizeData(object = rds, normalization.method = "LogNormalize",scale.factor = 10000)
rds <- FindVariableGenes(object = rds, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1,  y.cutoff = 1, do.contour=F)
use.gene = rds@var.genes
rds <- ScaleData(object = rds,vars.to.regress = c("nUMI", "percent.mito"),genes.use =use.gene)
rds <- RunPCA(object = rds, pc.genes = use.gene, do.print = FALSE)
rds <- FindClusters(object = rds, reduction.type = "pca", dims.use = 1:20, resolution = resolution, print.output = 0, save.SNN = TRUE,force.recalc = TRUE)
rds@meta.data[[res_str]] = as.numeric(rds@meta.data[[res_str]]) + 1
rds = SetAllIdent(rds, res_str)

# Run Non-linear dimensional reduction (tSNE)
rds <- RunTSNE(object = rds, dims.use = 1:20, do.fast = TRUE,check_duplicates = FALSE)
tryCatch({
  rds.markers <- FindAllMarkers(object = rds, genes.use = use.gene)
  rds.markers = dplyr::group_by(rds.markers,cluster) %>% dplyr::arrange(desc(avg_logFC))
}, error = function(e){
  print (paste0("no marker found: ", e))
  rds.markers <<- data.frame(cluster=double(),
                  gene=double(),
                  avg_logFC=double(),
                  pct.1=double(),
                  pct.2=double(),
                  p_val_adj=double())

})
print (rds.markers)
write_tsv(rds.markers,marker.out,col_names = T)

df.tsne = rds@dr$tsne@cell.embeddings
df.tsne = as.data.frame(df.tsne)
meta = rds@meta.data
dic = rds@meta.data[[res_str]]
names(dic) = rownames(rds@meta.data)
df.tsne$cluster = as.numeric(dic[rownames(df.tsne)])
df.gene = meta[,"nGene",drop=F]
colnames(df.gene) = "Gene_Counts"
df.all = cbind(df.tsne,df.gene)
write.table(df.all,tsne.out,sep="\t",col.names=NA,quote = F)

if (save_rds == 'True'){
  saveRDS(rds, rds.out)
}