# Introduction

This pipeline is currently mainly for TCR full-length assembly. The whole process includes 4 steps: sample, barcode, convert, assemble. The sample and barcode steps are the same as the previous process. For related information, please refer to the description of CeleScope. Here we mainly explain the steps of convert and assemble.

## Convert

### Features

- Convert Singleron barcode/UMI to 10X barcode/UMI.

### Input

- `--fq2` R2 reads from barcode step. Required.

### Output

- `{sample}/02.convert/barcode_correspond.txt` Recording barcodes correspondence. Two columns, first is Singleron barcode, second is 10X barcode.
- `{sample}/02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads with 10X barcode and UMI.
- `{sample}/02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads.

## Assemble

### Features

- TCR full length assembly.
- VDJ results match with scRNA-seq results.

### Input

- `--fqs_dir` Fastq file path after convert step. E.g. `02.convert`. Required.
- `--barcode_dic` File records barcodes correspondence. E.g. `{sample}/02.convert/barcode_correspond.txt`. Required.

### Other parameters

- `--species` Choose from ['hs', 'mmu']. Required.
- `--soft` Select cellranger version. ['3.0.2', '3.1.0', '4.0.0', '6.0.0']. Default `4.0.0`.
- `--mem` Memory for cellranger assembly. Default `10`.
- `--seqtype` Choose from ['TCR', 'BCR']. Required.
- `--match_dir` scRNA analysis directory. Required.

### Output

The assembly results are analyzed by cellranger and output in `{sample}/03.assemble/{sample}` directory. Replace the barcode in the cellranger with the Singleron barcode and output it to `{sample}/03.assemble/all` directory. Match with cell barcodes in the scRNA-seq, take the intersection, and output the barcode result to `{sample}/03.assemble/match` directory.

- `{sample}/03.assemble/all/filtered_contig_annotations.csv` 
- `{sample}/03.assemble/all/clonotypes.csv`
- `{sample}/03.assemble/all/filtered_contig.fasta`
- `{sample}/03.assemble/all/count.txt`
- `{sample}/03.assemble/match/match_contigs.csv`
- `{sample}/03.assemble/match/match_contig.fasta`
- `{sample}/03.assemble/match/match_clonotypes.csv`
- `{sample}/03.assemble/{sample}` Cellranger outs directory.
- `{sample}/03.assemble/{sample}_vdj_10X.sh` Cellranger run command.

# Example

Test path: 

`/SGRNJ03/randd/zhouxin/data/cellranger/`

Conda environment: 

`full_len_VDJ`

Run command: 

```
$ cat run.sh 
multi_vdj10X \
        --mapfile /SGRNJ03/randd/zhouxin/data/20210424.mapfile \
        --outdir ./ \
        --chemistry customized \
        --whitelist /SGRNJ02/RandD4/RD2019016/20200426/bclist_Reverse \
        --pattern U8C8L16C8L16C8L16 \
        --linker /SGRNJ02/RandD4/RD2019016/20200426/linker2 \
        --allowNoLinker --species mmu --thread 8 \
        --mod shell --soft 3.1.0 --seqtype TCR
```

```
$ cat shell/test1.sh 
set -eo pipefail
celescope vdj10X sample --outdir .//test1/00.sample --sample test1 --assay vdj10X --thread 8 --chemistry customized  --fq1 /SGRNJ03/randd/zhouxin/data/test_data/test1_1.fq 
celescope vdj10X barcode --outdir .//test1/01.barcode --sample test1 --assay vdj10X --thread 8 --chemistry customized --pattern U8C8L16C8L16C8L16 --whitelist /SGRNJ02/RandD4/RD2019016/20200426/bclist_Reverse --linker /SGRNJ02/RandD4/RD2019016/20200426/linker2 --lowNum 2 --allowNoLinker  --fq1 /SGRNJ03/randd/zhouxin/data/test_data/test1_1.fq --fq2 /SGRNJ03/randd/zhouxin/data/test_data/test1_2.fq 
celescope vdj10X convert --outdir .//test1/02.convert --sample test1 --assay vdj10X --thread 8 --method read2  --fq2 .//test1/01.barcode/test1_2.fq 
celescope vdj10X assemble --outdir .//test1/03.assemble --sample test1 --assay vdj10X --thread 8 --species mmu --soft 3.1.0 --mem 10 --seqtype TCR  --fqs_dir .//test1/02.convert --match_dir /SGRNJ03/randd/RD20040201_SCOPEv2_TCR/20210401_5/NJU_Pt_2_GOT/ --barcode_dic .//test1/02.convert/barcode_correspond.txt
```









