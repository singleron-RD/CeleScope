Both single-cell RNA-Seq and virus-enriched libraries were analyzed using CeleScope.

## Preprocessing
Barcodes and UMIs were extracted from R1 reads and corrected. Adapter sequences and poly A tails were trimmed from R2 reads using Cutadapt v3.7.

## Single-cell RNA-Seq library analysis
After proprocessing, R2 reads were aligned against the hg38 transcriptome using STAR v2.6.1a. Uniquely mapped reads were then assigned to exons with FeatureCounts(v2.0.1). Successfully Assigned Reads with the same cell barcode, UMI and gene were grouped together to generate the gene expression matrix for further analysis.

## Virus-enriched library analysis
After proprocessing, R2 reads were aligned against the EBV genome using STAR v2.6.1a with --outFilterMatchNmin set to 80. After obtaining the BAM file, in order to remove ambient contamination, we adopted a two-step filtration process. First, a UMI needs to be supported by a certain number of reads. Second, an EBV-positive cell needs to be supported by a certain number of UMIs. Thresholds for supporting reads and UMIs were determined by the Otsu's method.

## Reference
- MARTIN, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may 2011. ISSN 2226-6089.
- Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014
- Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21
- Bost, Pierre, et al. "Host-viral infection maps reveal signatures of severe COVID-19 patients." Cell 181.7 (2020): 1475-1488.
- Otsu, Nobuyuki. "A threshold selection method from gray-level histograms." IEEE transactions on systems, man, and cybernetics 9.1 (1979): 62-66.