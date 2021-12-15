# Introduction

This pipeline is currently mainly for TCR/BCR full-length assembly. The whole process includes 5 steps: sample, barcode, cutadapt, assemble, summarize. The sample, barcode, cutadapt steps are the same as the previous process. For related information, please refer to the description of CeleScope. Here we mainly explain the steps of assemble and summarize.

## Assemble

### Features

- Match with scRNA-seq reads.
- TCR/BCR full length assembly.

### Input

- `--fq2` R2 reads from cutadapt step. Required.
- `--match_dir' Match scRNA-seq directory. Required.

### Other parameters

- `--species` Choose from ['hg19', 'hg38', 'GRCm38']. Required.
- `--seqtype` Select TCR/BCR seq data type. Choose from ['TCR', 'BCR']. Required..
- `--barcodeRange` Barcode range in fq1, INT INT CHAR. Default `0 23 +`.
- `--umiRange` UMI range in fq1, INT INT CHAR. Default `24 -1 +`.
- `--trimLevel` INT: 0: no trim; 1: trim low quality; 2: trim unmatched. Default '1'.
- `--UMI_min` UMI number for cell determined, INT. Default 'auto'.

### Output

The assembly results are analyzed by Trust-4 and output in `{sample}/03.assemble/assemble` directory. Match with the scRNA-seq and output the result to `{sample}/03.assemble/match` directory.

- `03.assemble/match/{sample}_matched_R1.fq` New R1 reads matched with scRNA-seq
- `03.assemble/match/{sample}_matched_R1.fq` New R2 reads matched with scRNA-seq
- `03.assemble/assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.
- `03.assemble/assemble/{sample}_assign.out` read assignment results.
- `03.assemble/assemble/{sample}_assembled_reads.fa` Assembled raw reads.
- `03.assemble/assemble/{sample}_annot.fa` Assembled annotated contig sequences.
- `03.assemble/assemble/{sample}_full_len.fa` Assembled full length contig sequences.
- `03.assemble/assemble/report.out` Record assembled CDR3 types, read count and proportion of read count.
- `03.assemble/assemble/barcoderep.tsv` Record chain information for each barcode.
- `03.assemble/assemble/barcoderepfl.tsv` Record chain information for each barcode(preliminary filter).

## Summarize

### Features

- Filter cell number based on assembly result.
- Summarize result files and output.

### Input

- `--seqtype` Select TCR/BCR seq data type. Choose from ['TCR', 'BCR']. Required.
- `--min_read_count` Filter cell by read count number, int type required. Default '4'.

### Other parameters

- `--reads_assignment` File records reads assigned to contigs.. Required.
- `--assembled_fa` Reads used for assembly. Required..
- `--fq2` R2 reads from cutadapt step. Required.

### Output

- `04.summarize/clonetypes.tsv` High-level descriptions of each clonotype.
- `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
- `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
- `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig_annotations.csv.
- `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
- `04.summarize/{sample}_chain_filtered_contig.csv`Keep the 2 contigs with the highest UMI. This is a subset of filtered_contig.csv.
- `04.summarize/{sample}_chain_filtered_contig.fasta`Keep the 2 contigs with the highest UMI. This is a subset of filtered_contig.fasta.
- `04.summarize/{sample}_one_chain_contig.csv`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.csv.
- `04.summarize/{sample}_one_chain_contig.fasta`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.fasta.

# Example

Test path: 

`/SGRNJ03/randd/cjj/celedev/trust4/20211213_one_chain/3`

Conda environment: 

`vdj_trust4`

Run command: 

```
$ cat run.sh 
multi_trust_vdj \
        --mapfile ./test.mapfile \
        --outdir ./ \
        --chemistry customized \
        --whitelist /SGRNJ03/randd/cjj/zhouxin/zhouxin/data_base/bclist \
        --pattern U9C8L16C8L16C8L1 \
        --linker /SGRNJ03/randd/zhouxin/data_base/linker_4types_reversed_new \
        --allowNoLinker --species GRCm38 --thread 8 --mod sjm \
        --seqtype BCR \
```

```
$ cat shell/test.sh 
set -eo pipefail
celescope trust_vdj sample --outdir .//BJ_1116PzB_Auto_V2_Nlib/00.sample --sample BJ_1116PzB_Auto_V2_Nlib --assay trust_vdj --thread 8 --chemistry customized  --fq1 /SGRNJ03/DATA03/2112/20211202_4/R211126086_R1.fastq.gz
celescope trust_vdj barcode --outdir .//BJ_1116PzB_Auto_V2_Nlib/01.barcode --sample BJ_1116PzB_Auto_V2_Nlib --assay trust_vdj --thread 8 --chemistry customized --pattern U9C8L16C8L16C8L1 --whitelist /SGRNJ03/randd/cjj/zhouxin/zhouxin/data_base/bclist --linker /SGRNJ03/randd/zhouxin/data_base/linker_4types_reversed_new --lowNum 2 --allowNoLinker  --fq1 /SGRNJ03/DATA03/2112/20211202_4/R211126086_R1.fastq.gz --fq2 /SGRNJ03/DATA03/2112/20211202_4/R211126086_R2.fastq.gz
celescope trust_vdj cutadapt --outdir .//BJ_1116PzB_Auto_V2_Nlib/02.cutadapt --sample BJ_1116PzB_Auto_V2_Nlib --assay trust_vdj --thread 8 --minimum_length 20 --nextseq_trim 20 --overlap 10 --insert 150  --fq .//BJ_1116PzB_Auto_V2_Nlib/01.barcode/BJ_1116PzB_Auto_V2_Nlib_2.fq
celescope trust_vdj assemble --outdir .//BJ_1116PzB_Auto_V2_Nlib/03.assemble --sample BJ_1116PzB_Auto_V2_Nlib --assay trust_vdj --thread 8 --species GRCm38 --seqtype BCR --barcodeRange "0 23 +" --umiRange "24 -1 +" --trimLevel 1 --UMI_min auto  --fq2 .//BJ_1116PzB_Auto_V2_Nlib/02.cutadapt/BJ_1116PzB_Auto_V2_Nlib_clean_2.fq --match_dir /SGRNJ03/randd/RD20040201_SCOPEv2_TCR/20211205/BJ_1116PzB_Auto_V2_Nlib
celescope trust_vdj summarize --outdir .//BJ_1116PzB_Auto_V2_Nlib/04.summarize --sample BJ_1116PzB_Auto_V2_Nlib --assay trust_vdj --thread 8 --seqtype BCR --min_read_count auto  --reads_assignment .//BJ_1116PzB_Auto_V2_Nlib/03.assemble/assemble/BJ_1116PzB_Auto_V2_Nlib_assign.out --assembled_fa .//BJ_1116PzB_Auto_V2_Nlib/03.assemble/assemble/BJ_1116PzB_Auto_V2_Nlib_assembled_reads.fa --fq2 .//BJ_1116PzB_Auto_V2_Nlib/02.cutadapt/BJ_1116PzB_Auto_V2_Nlib_clean_2.fq
```









