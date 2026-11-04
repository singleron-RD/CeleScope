# scRNA-seq data processing (GEXSCOPE-V2 chemistry)

## Library structure and chemistry

Single-cell RNA-seq libraries were prepared with the Singleron GEXSCOPE-V2 chemistry. Cell barcodes in this chemistry are combinatorial: each barcode is assembled from three 9-nt segments (bc1, bc2 and bc3) drawn from three independent whitelists of 96, 96 and 192 sequences, respectively, separated by two 9-nt linkers, which yields a high-complexity barcode space for cell identification. Each transcript carries a 12-nt unique molecular identifier (UMI) and a poly(dT) capture sequence. The resulting read structure is:

| Read | Structure | Description |
|---|---|---|
| R1 | `C9-L9-C9-L9-C9-U12` | three cell-barcode segments, two linkers, 12-nt UMI |
| R2 | cDNA insert | transcript sequence, poly(A) read-through removed by 3′ adapter clipping |

## Read processing, alignment and quantification

Raw sequencing reads were demultiplexed and processed with CeleScope v2.15.4 (`multi_rna`; https://github.com/singleron-RD/CeleScope). Reads were first trimmed of the 3′ poly(A) adapter, and were aligned to the Ensembl release 110 human reference genome (GRCh38) with STAR (v2.7.11a) in STARsolo mode (`https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md`). The reference index was built with `celescope rna mkref` from the unfiltered GTF annotation so that short non-coding RNAs were retained.

Cell barcodes and UMIs were extracted in STARsolo using the complex barcode layout (`--soloType CB_UMI_Complex`), with the three barcode segments and the UMI specified at their read positions; the two 9-nt linkers were supplied as the expected adapter sequence to validate correct library structure. Barcodes were matched against the predefined whitelists with `--soloCBmatchWLtype EditDist_2`, and UMIs were collapsed to remove PCR duplicates. Only alignments of at least 50 matched bases were retained (`--outFilterMatchNmin 50`).

## Cell calling

Empty droplets were removed with STARsolo cell calling using the EmptyDrops implementation of CellRanger (`--soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000`), which distinguishes real cells from ambient RNA background on the basis of the full barcode count distribution; barcodes passing the filter were retained as cells. The same run produced raw and filtered (cell-associated) barcode × gene matrices.

## Gene quantification

Counts were generated simultaneously at gene level and with the pre-mRNA-permissive `GeneFull_Ex50pAS` strategy (`--soloFeatures GeneFull_Ex50pAS Gene`), in which a read is assigned to a gene when at least 50% of the read overlaps an exon, thereby retaining intronic reads from unspliced transcripts. The `GeneFull_Ex50pAS` counts were used for the HTML report and for all downstream analyses, as they improve sensitivity for single-nucleus and low-complexity libraries.

## Downstream analysis

The filtered count matrix was analysed in Scanpy v1.9.3 (AnnData v0.7.6). Per-cell quality-control metrics were calculated, including the number of detected genes, total UMI counts and the fraction of counts from mitochondrial genes (genes with the `MT-` prefix), which were reported per library without hard filtering. Counts were library-size normalized to 10,000 counts per cell, log1p-transformed and stored as a `normalized` layer that was used for all subsequent steps and for marker detection.

Highly variable genes were selected with the `seurat` flavor (`min_disp = 0.5`, `min_mean = 0.0125`, `max_mean = 3`, 20 bins), the data were scaled with zero centering and a maximum value of 10, and principal component analysis was run with 50 components on the highly variable genes. A 15-nearest-neighbor graph was constructed on the top 25 principal components and used for both t-SNE and UMAP embedding (UMAP with `min_dist = 0.5`).

Cells were clustered with the Leiden algorithm at resolution 1.2. Cluster marker genes were identified by Wilcoxon rank-sum tests (`sc.tl.rank_genes_groups`, one cluster versus all remaining cells) on the normalized layer, with Benjamini–Hochberg correction for multiple testing. Genes with adjusted *P* < 0.05 and positive log2 fold change were reported, retaining the 100 most significant markers per cluster. The analysis output was exported as an `.h5ad` object retaining raw and normalized expression, cluster assignments, t-SNE and UMAP coordinates, and QC metrics. Where requested, cell-type annotation was additionally performed with CellTypist.

## Reference

- Kaminow B, Yunusov D, Dobin A. STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data. bioRxiv. 2021.05.05.442755; doi: https://doi.org/10.1101/2021.05.05.442755

- Dyer SC, Austine-Orimoloye O, Azov AG, Barba M, Barnes I, Barrera-Enriquez VP, et al. Ensembl 2025. Nucleic Acids Research. 2025;53(D1):D948–D957. doi: https://doi.org/10.1093/nar/gkae1071

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19(1):15. doi: https://doi.org/10.1186/s13059-017-1382-0