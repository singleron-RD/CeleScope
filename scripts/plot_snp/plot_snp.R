library(tidyverse)
library(argparser)

# --plot setting ----------------------------------------------------------------------

DEVICE = "pdf"
WIDTH = 7
HEIGHT = 7

# --function ----------------------------------------------------------------------

tsne_add_genetype <- function(variant, gt_df, tsne_df){
  type <- gt_df %>% 
    filter(str_detect(pos,str_trim(variant))) %>% 
    rbind(colnames(.)) %>% 
    select(-pos) %>% 
    t() %>% 
    as_tibble(.name_repair = 'unique')
   names(type) <- c("genotypes", "barcode")
   tsne_gene_type <- tsne_df %>% left_join(type, by = "barcode")
   return(tsne_gene_type)
}

protein_to_gtpos <- function(protein, variant_df){
    variant_gt <- c()
    variant_select <- variant_df %>% 
        filter(str_detect(Protein,str_trim(protein)))
    chrom <- variant_select$Chrom
    pos <- variant_select$Pos
    id <- variant_select$id
    variant <- unique(str_c(chrom, "_", pos, "_", id))
    return(variant)
}

position_to_gtpos <- function(position, gt_df){
    v_gt <- gt_df %>% filter(str_detect(pos, position)) %>% .$pos
    return(v_gt)
}

plot_gttsne <- function(variant, gt_df, tsne_df){
    df_tsne <- tsne_add_genetype(variant, gt_df, tsne_df)
    p <- ggplot(df_tsne, mapping = aes(tSNE_1,tSNE_2, color = genotypes)) + 
                geom_point()
    return(p)
}

plot_position <- function(variant_input, sample, gt_df, tsne_df){
    for (variant in variant_input) {
        v_gt <- position_to_gtpos(variant, gt_df)
        for (v in v_gt){
            p <- plot_gttsne(v, gt_df, tsne_df)
            p <- p + ggtitle(str_glue("position: {v}"))
            ggsave(filename = str_glue("{out_dir}/{sample}_{v}.{DEVICE}"),plot = p,device = DEVICE,width = WIDTH,height = HEIGHT)
        }
    }
}

plot_protein <- function(protein_input, sample, gt_df, variant_df, tsne_df){
    for (protein in protein_input){
        v_gt <- protein_to_gtpos(protein, variant_df)
        for (v in v_gt){
            p <- plot_gttsne(v, gt_df, tsne_df)
            p <- p + ggtitle(str_glue("protein / position: {protein} / {v}"))
            ggsave(filename = str_glue("{out_dir}/{sample}_{protein}_{v}.{DEVICE}"),plot = p,device = DEVICE,width = WIDTH,height = HEIGHT)
        }
    }
}

# --argparser ----------------------------------------------------------------------

argv <- arg_parser('')
argv <- add_argument(argv,"--snp_dir", short = "-d", 
                     help="CeleScope snp directory")
argv <- add_argument(argv,"--out_dir", short = "-o", 
                     help="Output directory.", default = "./")
argv <- add_argument(argv,"--protein", default = "", 
                     help="Optional. Plot variants with these protein changes. Multiple protein changes are separated by comma.")
argv <- add_argument(argv,"--position", default = "", 
                    ,help="Optional. Plot variants at these positions. Multiple positions are separated by comma.")
argv <- parse_args(argv)

# -- load data ----------------------------------------------------------------------

snp_dir <- argv$snp_dir
out_dir <- argv$out_dir

if(dir.exists(out_dir) == FALSE) {
  dir.create(out_dir)
}

sample <- str_split(snp_dir, "/")[[1]]
if (sample[length(sample)] == ""){
    sample <- sample[length(sample) - 1]
}else{
    sample <- sample[length(sample)]
}

tsne <- str_glue("{snp_dir}/09.analysis_snp/{sample}_tsne_coord.tsv")
gt <- str_glue("{snp_dir}/09.analysis_snp/{sample}_gt.csv")
variant <- str_glue("{snp_dir}/09.analysis_snp/{sample}_variant_table.csv")

protein <- str_split(argv$protein, ",")[[1]]
position <- str_split(argv$position, ",")[[1]]

gt_df <- read_csv(gt) %>% rename("pos" = "X1")
variant_df <- read_csv(variant) %>% mutate(id = c(1:nrow(.)))
tsne_df <- read_tsv(tsne) %>% rename("barcode" = "X1")

# -- plot ----------------------------------------------------------------------

if (any(position == "") & any(protein != "")){
    plot_protein(protein, sample, gt_df, variant_df, tsne_df)
}else if(any(position != "") & any(protein == "")){
    plot_position(position, sample, gt_df, tsne_df)
}else if(any(position != "") & any(protein != "")){
    plot_position(position, sample, gt_df, tsne_df)
    plot_protein(protein, sample, gt_df, variant_df, tsne_df)
}else{
    message("None input of protein or position!")
}