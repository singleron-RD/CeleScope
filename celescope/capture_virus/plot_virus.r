library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--tsne_capture_virus", help="tsne_capture_virus file")
argv <- add_argument(argv,"--outdir", help="output dir")
argv <- parse_args(argv)

tsne_capture_virus = argv$tsne_capture_virus
outdir = argv$outdir

df = read_tsv(tsne_capture_virus)
ncol = dim(df)[2]
cols = colnames(df[,c(6:ncol)])

for (col in cols){
    #df[col] = apply(df[col],FUN=str,MARGIN=1)
    p.fea = ggplot(df,aes_string(x="tSNE_1",y="tSNE_2",color=col)) + 
        geom_point() + 
        scale_color_gradient(low = "gray", high = "blue")
    p.vio = ggplot(df,aes_string(x="as_factor(cluster)",y=col,fill="as_factor(cluster)")) + 
        geom_violin() + labs(x = "cluster", y = col, group = "cluster")
    
    pdf.name.fea = paste0(outdir,"/",col,"_capture_virus_featureplot.pdf")
    pdf.name.vio = paste0(outdir,"/",col,"_capture_virus_violinplot.pdf")
    pdf(pdf.name.fea,height=10,width=12)
    print (p.fea)

    pdf(pdf.name.vio,height=10,width=12)
    print (p.vio)

    dev.off()
}
