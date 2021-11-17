library(vcfR)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv, "--vcf", help="vcf file path")
argv <- add_argument(argv, "--out", help="output file path")
argv <- add_argument(argv, "--sample", help="sample name")

argv <- parse_args(argv)

vcf = read.vcfR(argv$vcf, verbose = FALSE)
gt <- extract.gt(vcf)

write.csv(gt, argv$out, quote=FALSE)


