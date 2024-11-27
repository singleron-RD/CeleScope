## Preprocessing
Barcodes and UMIs were extracted from R1 reads and corrected. Adapter sequences and poly A tails were trimmed from R2 reads using Cutadapt v3.7.

## Virus-enriched library analysis
After preprocessing, R2 reads were aligned against the EBV genome(https://www.ncbi.nlm.nih.gov/nuccore/NC_007605) using STAR v2.6.1a with --outFilterMatchNmin set to 80. Uniquely mapped reads were then assigned to EBV genes with FeatureCounts(v2.0.1). After obtaining the BAM file, to remove ambient contamination, we adopted a filtering step borrowed from Viral-Track(https://github.com/PierreBSC/Viral-Track). A valid UMI needs to be supported by a certain number of reads. The threshold for supporting reads was determined by the Otsu's method.

## Reference
- MARTIN, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may 2011. ISSN 2226-6089.
- Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014
- Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21
- Bost, Pierre, et al. "Host-viral infection maps reveal signatures of severe COVID-19 patients." Cell 181.7 (2020): 1475-1488.
- Otsu, Nobuyuki. "A threshold selection method from gray-level histograms." IEEE transactions on systems, man, and cybernetics 9.1 (1979): 62-66.