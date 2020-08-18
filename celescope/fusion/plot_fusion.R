library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--tsne_fusion", help="tsne_fusion file")
argv <- add_argument(argv,"--outdir", help="output dir")
argv <- parse_args(argv)

tsne_fusion = argv$tsne_fusion
outdir = argv$outdir

df = read_tsv(tsne_fusion)
ncol = dim(df)[2]
cols = colnames(df[,c(6:ncol)])
str = function(x){ifelse(x>0,"Fusion","No Fusion")}
for (col in cols){
    df[col] = apply(df[col],FUN=str,MARGIN=1)
    p = ggplot(df,aes_string(x="tSNE_1",y="tSNE_2",color=col)) + geom_point()
    pdf.name = paste0(outdir,"/",col,"_fusion.pdf")
    pdf(pdf.name,height=10,width=12)
    print (p)
    dev.off()
}
