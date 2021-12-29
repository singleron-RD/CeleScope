library(DropletUtils)
library(Seurat)
library(argparser)
library(tidyverse)

argv <- arg_parser('')
argv <- add_argument(argv,"--outdir", help="the output dir.",default=getwd())
argv <- add_argument(argv,"--sample", help="sample")
argv <- add_argument(argv,"--matrix_dir", help="matrix_dir")
argv <- add_argument(argv,"--threshold", help="UMI threshold")
argv <- parse_args(argv)

#read args
outdir = argv$outdir
sample = argv$sample
matrix_dir = argv$matrix_dir
threshold = as.numeric(argv$threshold)

mtx = Seurat::Read10X(matrix_dir)
bcrank = barcodeRanks(mtx, lower=threshold/2)
inflection = as.numeric(bcrank@metadata$inflection)
knee = as.numeric(bcrank@metadata$knee)

origin.cell = sum(bcrank$total >= threshold)
inflection.cell = sum(bcrank$total >= inflection)
rescued.cell = inflection.cell - origin.cell
if (rescued.cell < 0){
    rescued.cell = 0
}
print(paste('rescued cell number:',rescued.cell))

df = tibble(
    inflection=as.character(inflection),
    knee=as.character(knee),
    rescued_cell=as.character(rescued.cell)
)
print(df)
df.out = stringr::str_glue('{outdir}/{sample}_rescue.tsv')
write_tsv(df, df.out)





