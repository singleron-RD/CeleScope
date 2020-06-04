library(Seurat)
rds = Read10X(data.dir = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/saturation/run_low/TACTTCGG/mtx", 
              gene.column = 1)

rds = readRDS("/SGRNJ02/RandD4/RD2020001_SCOPEV2/20200528/result4/rds/all_TSNE_samples.rds")
head(rds@meta.data)

##

df.tsne = rds@dr$tsne@cell.embeddings
df.tsne = as.data.frame(df.tsne)
dic = rds@meta.data$res.0.6
names(dic) = rownames(rds@meta.data)
df.tsne$cluster = as.character(dic[rownames(df.tsne)])
write.table(df.tsne,"/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/scope_tools_1.0/out/06.analysis/tsne_cluster.tsv",
            sep="\t",col.names=NA,quote = F)
