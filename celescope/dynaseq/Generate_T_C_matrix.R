args <- commandArgs(T)

require("reshape2")
require("tidyr")
require("dplyr")
require("Matrix")

my.count1 <- read.table(args[1],h=F)

my.count1$V1 <- as.character(my.count1$V1)
my.count1$gene <- gsub("--C","",my.count1$V1)
my.count1$gene <- gsub("--T","",my.count1$gene)
cells.keep <- my.count1 %>% dplyr::distinct(V2,V3,gene) %>% group_by(V2) %>% dplyr::summarize(count=n()) %>% arrange(desc(count)) %>% .$V2 %>% as.character 

inds <- as.numeric(args[2])
if (length(cells.keep) > inds) {
 cells.keep2 <- head(cells.keep,inds)
}else{ cells.keep2 <- cells.keep}

my.count1 <- my.count1 %>% filter(V2 %in% cells.keep2) %>% droplevels
my.count1$type <- "C"
my.count1[grep("--T",my.count1$V1),]$type <- "T"
my.count2 <- dcast(my.count1,gene+V2+V3 ~ type, value.var = "V4")
my.count2[is.na(my.count2)] <- 0

if(! "C" %in% colnames(my.count2))
{
  my.count2$C <- 0;
}
if (ncol(my.count2) !=5) {
    stop("Error! Please verify the count data frame!\n");
}
my.count2 <- my.count2 %>% arrange(gene,V2,V3,C,T)
my.count2 %>% mutate(type = ifelse(C > 0,"C","T")) -> my.count2 
my.count3 <- my.count2 %>% group_by(gene,type,V2) %>% dplyr::summarize(count=n())
my.count3$gene2 <- paste(my.count3$gene,my.count3$type,sep="--")
my.count3$V2 <- as.factor(my.count3$V2)
my.count3$gene2 <- as.factor(my.count3$gene2)
data.sparse = sparseMatrix(as.integer(my.count3$gene2), as.integer(my.count3$V2), x = my.count3$count)
colnames(data.sparse) = levels(my.count3$V2)
rownames(data.sparse) = levels(my.count3$gene2)
ord <- sort(colSums(data.sparse),decreasing = T)
data.sparse <- data.sparse[,names(ord)]
saveRDS(data.sparse,file=args[3])
outtsv<-paste(args[3],"tsv", sep = ".")
write.table(as.matrix(data.sparse), file = outtsv, sep = "\t", quote = F, row.names = T)





