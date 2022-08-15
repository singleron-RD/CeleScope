Raw reads were processed using CeleScope v1.5.1(https://github.com/singleron-RD/CeleScope/tree/v1.5.1). Briefly, Barcodes and UMIs were extracted from R1 reads and corrected. Adapter sequences and poly A tails were trimmed from R2 reads and the trimmed R2 reads were aligned against the GRCh38 transcriptome using STAR v2.6.1b. Then, Bcftools is used for variant calling, with default parameters. We borrowed the method in the Viral-Track to filter the background contamination. For each variant, we used the Otsu's method to determine the supporting reads threshold T. If the number of supporting reads for a variant in a cell is less than T, it was set to zero.

## Reference
- Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21
- Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93.
- Bost, Pierre, et al. "Host-viral infection maps reveal signatures of severe COVID-19 patients." Cell 181.7 (2020): 1475-1488.
- Otsu, Nobuyuki. "A threshold selection method from gray-level histograms." IEEE transactions on systems, man, and cybernetics 9.1 (1979): 62-66.
