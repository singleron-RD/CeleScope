library(Seurat)
library(tidyverse)
rds = Read10X(data.dir = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/saturation/run_low/TACTTCGG/mtx", 
              gene.column = 1)

rds = readRDS("/SGRNJ02/RandD4/RD2020001_SCOPEV2/20200528/result4/rds/all_TSNE_samples.rds")
head(rds@meta.data)
meta = rds@meta.data
##

df.tsne = rds@dr$tsne@cell.embeddings
df.tsne = as.data.frame(df.tsne)
dic = rds@meta.data$res.0.6
names(dic) = rownames(rds@meta.data)
df.tsne$cluster = as.character(dic[rownames(df.tsne)])

df.gene = meta[,"nGene",drop=F]
colnames(df.gene) = "Gene_Counts"
df.all = cbind(df.tsne,df.gene)

#rds = RunUMAP(rds)

write.table(df.all,"/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/scope_tools_1.0/out/06.analysis/tsne_coord.tsv",
            sep="\t",col.names=NA,quote = F)



##markers
library(tidyverse)
markers.df = read_csv("/SGRNJ02/RandD4/RD2020001_SCOPEV2/20200528/result4/csv/all_markers.csv")
markers.df$cluster = markers.df$cluster+1
markers.df
top.df = markers.df %>% group_by(cluster) %>% top_n(20,avg_logFC)
write_tsv(top.df,"/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/scope_tools_1.0/out/06.analysis/markers.tsv")
