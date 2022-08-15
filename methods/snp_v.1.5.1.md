Both single-cell RNA-Seq and panel-enriched libraries were analyzed using CeleScope.

## Preprocessing
Barcodes and UMIs were extracted from R1 reads and corrected. Adapter sequences and poly A tails were trimmed from R2 reads using Cutadapt v3.7. 

## Single-cell RNA-Seq library analysis
After proprocessing, R2 reads were aligned against the hg38 transcriptome using STAR v2.6.1b. Uniquely mapped reads were then assigned to exons with FeatureCounts(v2.0.1). Successfully Assigned Reads with the same cell barcode, UMI and gene were grouped together to generate the gene expression matrix for further analysis.

## Single-cell panel-enriched library analysis
After proprocessing, R2 reads were aligned against the hg38 transcriptome using STAR v2.6.1b and only reads with valid cell barcode were kept. Then, Bcftools is used for variant calling, with default parameters. We borrowed the method in the Viral-Track to filter the background contamination. For each variant, we used the Otsu's method to determine the supporting reads threshold T. If the number of supporting reads for a variant in a cell is less than T, it was set to zero.

## Reference
-- MARTIN, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may 2011. ISSN 2226-6089.
- Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21
- Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93.
- Bost, Pierre, et al. "Host-viral infection maps reveal signatures of severe COVID-19 patients." Cell 181.7 (2020): 1475-1488.
- Otsu, Nobuyuki. "A threshold selection method from gray-level histograms." IEEE transactions on systems, man, and cybernetics 9.1 (1979): 62-66.
