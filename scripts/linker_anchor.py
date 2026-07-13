import argparse

import pysam

from celescope.tools import utils


LINKER1 = "ACGATG"
LINKER2 = "CATAGT"
BC_LEN = 9
LINKER_LEN = 6


def find_linker_pos(seq):
    """
    Find linker position in seq by searching linker sequences directly.
    Pattern: C9L6C9L6C9U12
    Linkers: ACGATG, CATAGT
    No mismatch allowed for linkers.
    Returns (first_linker_start, second_linker_start) if both found, else (-1, -1).
    """
    first_linker_start = seq.find(LINKER1)
    if first_linker_start == -1:
        return -1, -1

    second_linker_start = seq.find(LINKER2, first_linker_start + LINKER_LEN)
    if second_linker_start == -1:
        return -1, -1

    return first_linker_start, second_linker_start


def extract_bc_umi(seq, first_linker_start, second_linker_start):
    """
    Extract barcode and UMI from seq given linker positions.
    Pattern: C9L6C9L6C9U12
    Returns (barcode, umi) if bc2 length is valid, else (None, None)
    """
    bc2_start = first_linker_start + LINKER_LEN
    bc2 = seq[bc2_start:second_linker_start]
    if len(bc2) != BC_LEN:
        return None, None

    bc1_start = first_linker_start - BC_LEN
    bc1 = seq[bc1_start:first_linker_start]

    bc3_start = second_linker_start + LINKER_LEN
    bc3 = seq[bc3_start : bc3_start + BC_LEN]

    umi_start = bc3_start + BC_LEN
    umi = seq[umi_start : umi_start + 12]

    barcode = f"{bc1}_{bc2}_{bc3}"
    return barcode, umi


def run(fq1, fq2, out_fq2):
    total_reads = 0
    valid_reads = 0

    with pysam.FastxFile(fq1) as fh1, pysam.FastxFile(fq2) as fh2:
        with utils.generic_open(out_fq2, "wt") as out_fh:
            for read1, read2 in zip(fh1, fh2):
                total_reads += 1
                seq1 = read1.sequence

                first_linker_start, second_linker_start = find_linker_pos(seq1)
                if first_linker_start == -1:
                    continue

                barcode, umi = extract_bc_umi(
                    seq1, first_linker_start, second_linker_start
                )
                if barcode is None:
                    continue
                valid_reads += 1

                read_name = f"{barcode}:{umi}:{total_reads}"
                out_fh.write(utils.fastq_line(read_name, read2.sequence, read2.quality))

    print(f"Total reads: {total_reads}")
    print(f"Valid reads: {valid_reads}")
    print(f"Valid rate: {valid_reads / total_reads * 100:.2f}%")


def main():
    parser = argparse.ArgumentParser(
        description="Extract barcode and UMI from R1 using linker anchor, add to R2 read name."
    )
    parser.add_argument("--fq1", help="R1 fastq file", required=True)
    parser.add_argument("--fq2", help="R2 fastq file", required=True)
    parser.add_argument("--out_fq2", help="Output R2 fastq file", required=True)
    args = parser.parse_args()
    run(args.fq1, args.fq2, args.out_fq2)


if __name__ == "__main__":
    main()
