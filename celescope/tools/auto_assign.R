library(argparser)
library(Seurat)
library(tidyverse)

argv <- arg_parser('')
argv <- add_argument(argv,"--outdir", help="the output dir.",default=getwd())
argv <- add_argument(argv,"--sample", help="sample")
argv <- add_argument(argv,"--rds",help="Seurat rds")
argv <- add_argument(argv,"--type_marker_tsv",help="cell type marker tsv")
#argv <- add_argument(argv,"--resolution", help="tSNE resolution",default=0.8)
argv <- parse_args(argv)

#read args

outdir <- argv$outdir
sample <- argv$sample
rds <- argv$rds
type_marker_tsv <- argv$type_marker_tsv
#resolution <- argv$resolution
rds <- argv$rds
#origin.cluster <- paste("res.",resolution,sep="")

print ("reading RDS")
all_data <- readRDS(rds)
print ("done.")
marker_file <- read_tsv(type_marker_tsv)

#
cell_name <- (marker_file)[,1,drop=T]
n_cell_name <- length(cell_name)

#reset
#all_data <- SetAllIdent(object = all_data, id = origin.cluster)
clusters <- sort(unique(all_data@active.ident))

#create dir
auto_dir <- stringr::str_glue('{outdir}/{sample}_auto_assign/')
png_dir <- stringr::str_glue('{auto_dir}/{sample}_png/')
dir.create(auto_dir)
dir.create(png_dir)

# type marker
c = 0 
for (cluster in clusters){
  index = 0
  for (cell in cell_name){
    index = index + 1
    pos = unlist(strsplit(marker_file[index,2,drop=T],","))
    neg = tryCatch(unlist(strsplit(marker_file[index,3,drop=T],",")) ,error=function(e){} )
    for (F in pos){
      tryCatch({
        dat <- FindMarkers(all_data,feature=F,ident.1=cluster,min.pct = 0,logfc.threshold = -Inf)
        dat$cell_type <- cell
        dat$cluster <- cluster
        dat <- rownames_to_column(dat,var="gene")
        dat$type <- "positive"
        if (c==0){
          all_dat <- dat
          c = c + 1
        } else {
          all_dat <- rbind(all_dat,dat)
          }
        }
        ,error=function(e){print(paste0(F," not found in cluster ",cluster)) })
    }

    if (!is.na(neg) && !is.null(neg)){
    	for (F in neg){
      	tryCatch({
        dat <- FindMarkers(all_data,feature=F,ident.1=cluster,min.pct = 0,logfc.threshold = -Inf)
        dat$cell_type <- cell
        dat$cluster <- cluster
        dat <- rownames_to_column(dat,var="gene")
        dat$type <- "negative"
        if (c==0){
          all_dat <- dat
          c = c + 1
        } else {
          all_dat <- rbind(all_dat,dat)
          }
        }
        ,error=function(e){print(paste0(F," not found in cluster ",cluster)) })
    	}
    }

  }
}

all_dat$cluster <- as.numeric(all_dat$cluster) + 1
all_dat <- mutate(all_dat,pct.diff=pct.1-pct.2)
exp.out = stringr::str_glue('{auto_dir}/{sample}_type_marker_exp.tsv')
write_tsv(all_dat, exp.out)

# plot
color2 <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple",
"DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4",
"#4682B4","#FFDAB9","#708090","#836FFF","#CDC673",
"#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80",
"#6A5ACD","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62",
"#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9","Green")

exp <- all_dat
a <- group_by(exp,cluster,cell_type)
for (cluster in clusters){
  c = a[a$cluster==cluster,]

  png(paste0(png_dir,cluster,"_pctdiff.png"),width=1200,height=1000)
  p1 <- ggplot(c,aes(x=interaction(gene,cell_type,type),y=pct.diff,fill=cell_type)) +geom_bar(stat="identity")+ coord_flip() + scale_fill_manual(values=color2)
  print (p1)
  dev.off()

  png(paste0(png_dir,cluster,"_logfc.png"),width=1200,height=1000)
  p2 <- ggplot(c,aes(x=interaction(gene,cell_type,type),avg_log2FC,fill=cell_type)) +geom_bar(stat="identity")+ coord_flip() + scale_fill_manual(values=color2)
  print (p2)
  dev.off()
}

# auto assign
exp[exp$type=="negative",]$avg_log2FC = -(exp[exp$type=="negative",]$avg_log2FC)
exp[exp$type=="negative",]$pct.diff = -(exp[exp$type=="negative",]$pct.diff)
a <- group_by(exp,cluster,cell_type)
as <- summarize(a,avg_pct.diff=mean(pct.diff),avg_log2FC=mean(avg_log2FC),max_p_val_adj=max(p_val_adj))
as1 <- group_by(ungroup(as),cluster)
as1 <- mutate(as1,pct_rank = rank(avg_pct.diff),
              logfc_rank= rank(avg_log2FC),total_rank=pct_rank+logfc_rank)
as2 <- as1 %>% ungroup %>% group_by(cluster) %>% 
  filter(total_rank==max(total_rank)) %>% arrange(as.numeric(cluster))
as3 <- select(as2,cluster,cell_type,avg_pct.diff,avg_log2FC,max_p_val_adj)
as3[(as3$avg_pct.diff < 0 | as3$avg_log2FC < 0),]$cell_type = 'NA'
res.out = stringr::str_glue('{auto_dir}/{sample}_auto_cluster_type.tsv')
write_tsv(as3, res.out)