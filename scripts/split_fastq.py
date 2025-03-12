import os
import gzip
import argparse
import pysam
from collections import defaultdict


def load_barcodes(barcode_file):
    barcodes = set()
    with gzip.open(barcode_file, "rt") as f:
        for line in f:
            barcodes.add(line.strip())  # 移除换行符并添加到集合
    return barcodes


def split_fastq_by_bc(input_fastq, barcode_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)  # 确保输出目录存在
    valid_barcodes = load_barcodes(barcode_file)  # 读取 barcodes
    print(f"Loaded {len(valid_barcodes)} barcodes from {barcode_file}")

    barcode_reads = defaultdict(list)  # 存储 barcode 对应的 reads

    # 读取 FASTQ 文件
    with pysam.FastxFile(input_fastq) as fq:
        for entry in fq:
            name_parts = entry.name.split(":")  # 拆分 read name
            if len(name_parts) >= 1:
                bc = name_parts[0]  # 提取 barcode
                if bc in valid_barcodes:  # 只处理指定的 barcodes
                    fastq_record = (
                        f"@{entry.name}\n{entry.sequence}\n+\n{entry.quality}\n"
                    )
                    barcode_reads[bc].append(fastq_record)

    # 按 barcode 写入 gzip 压缩的 FASTQ 文件
    for bc, reads in barcode_reads.items():
        output_path = os.path.join(output_dir, f"{bc}.fastq.gz")
        with gzip.open(output_path, "wt") as out_fq:
            out_fq.writelines(reads)

    print(f"Splitting completed! {len(barcode_reads)} matched barcodes found.")


def main():
    parser = argparse.ArgumentParser(
        description="Split FASTQ file based on barcode (BC) in read names"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input FASTQ file (can be gzipped)"
    )
    parser.add_argument(
        "-b",
        "--barcode",
        required=True,
        help="Barcode file (gzip compressed, one barcode per line)",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output directory for split FASTQ files"
    )

    args = parser.parse_args()

    split_fastq_by_bc(args.input, args.barcode, args.output)


if __name__ == "__main__":
    main()
