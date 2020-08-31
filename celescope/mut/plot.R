library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--tsne_mut", help="tsne_mut file")
argv <- add_argument(argv,"--outdir", help="output dir")
argv <- parse_args(argv)

tsne_mut = argv$tsne_mut
outdir = argv$outdir

df = read_tsv(tsne_mut)
ncol = dim(df)[2]
cols = colnames(df[,c(6:ncol)])
str = function(x){ifelse(x>0,"Mutation","WT")}
for (col in cols){
    df[col] = apply(df[col],FUN=str,MARGIN=1)
    p = ggplot(df,aes_string(x="tSNE_1",y="tSNE_2",color=col)) + geom_point() + labs(title=col)
    pdf.name = paste0(outdir,"/",col,"_mut.pdf")
    pdf(pdf.name,height=10,width=12)
    print (p)
    dev.off()
}
