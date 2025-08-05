"""
Split sample fastqs based on sample barcode on R2 start position
"""

import pysam
import argparse
import gzip
import os

from celescope.tools import utils


def get_barcode_sample_dict(barcode_sample_file):
    barcode_sample_dict = utils.two_col_to_dict(barcode_sample_file)
    bcs = list(barcode_sample_dict.keys())
    for bc in bcs:
        if len(bc) != len(bcs[0]):
            raise ValueError("All barcodes must have the same length")
    return barcode_sample_dict, len(bcs[0])


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fq1", required=True, help="Input fastq file 1 (R1)")
    parser.add_argument("--fq2", required=True, help="Input fastq file 2 (R2)")
    parser.add_argument(
        "-b", "--barcode_sample", required=True, help="barcode sample tsv file."
    )
    parser.add_argument("-o", "--output", help="Output folder", default="./outs")
    return parser.parse_args()


def main():
    args = get_args()
    barcode_sample, barcode_length = get_barcode_sample_dict(args.barcode_sample)
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    sample_fh = {}

    fq1 = pysam.FastxFile(args.fq1)
    fq2 = pysam.FastxFile(args.fq2)
    for read1, read2 in zip(fq1, fq2):
        bc = read2.sequence[:barcode_length]
        sample = barcode_sample.get(bc, "other")
        if sample not in sample_fh:
            fq1_out = gzip.open(
                f"{output_dir}/{sample}_R1.fastq.gz", "wt", compresslevel=1
            )
            fq2_out = gzip.open(
                f"{output_dir}/{sample}_R2.fastq.gz", "wt", compresslevel=1
            )
            sample_fh[sample] = (fq1_out, fq2_out)
        sample_fh[sample][0].write(str(read1) + "\n")
        sample_fh[sample][1].write(str(read2) + "\n")
    fq1.close()
    fq2.close()


if __name__ == "__main__":
    main()
