Spatial data were processed using the CeleScope pipeline (v2.11.0; https://github.com/singleron-RD/CeleScope). Specifically, raw reads were demultiplexed and aligned to the Ensembl release 110 human reference genome via STARsolo (v2.7.11b) to generate the gene expression matrix. Histological images were cropped and registered to spatial barcodes using AtlasXbrowser(https://github.com/singleron-RD/AtlasXbrowser), generating the final spot-level coordinates for tissue mapping.

## Reference

- Kaminow B, Yunusov D, Dobin A. STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data. bioRxiv. 2021.05.05.442755; doi: https://doi.org/10.1101/2021.05.05.442755

- Dyer SC, Austine-Orimoloye O, Azov AG, Barba M, Barnes I, Barrera-Enriquez VP, et al. Ensembl 2025. Nucleic Acids Research. 2025;53(D1):D948–D957. doi: https://doi.org/10.1093/nar/gkae1071