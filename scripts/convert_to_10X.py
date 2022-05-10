import pandas as pd
import pysam
import os
import argparse
from xopen import xopen

UMI_10X_LEN = 10
TSO = "TTTCTTATATGGG"
SEQ_LEN = 150


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")


def get_10X_barcode(cellranger_path):
    BARCODE_10X_FILE = os.path.dirname(cellranger_path) + "/lib/python/cellranger/barcodes/737K-august-2016.txt"
    return BARCODE_10X_FILE


def fastq_line(name, seq, qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'


def convert_seq(sgr_barcode, umi, barcode_dict, barcodes_10X, seq2, qual2):
    '''
    barcode_dict - key:SGR barcode; value:10X barcode
    '''
    if sgr_barcode in barcode_dict:
        barcode_10X = barcode_dict[sgr_barcode]
    else:
        barcode_10X = barcodes_10X.readline().strip()
        barcode_dict[sgr_barcode] = barcode_10X

    if len(umi) > UMI_10X_LEN:
        umi_10X = umi[0:UMI_10X_LEN]
    elif len(umi) < UMI_10X_LEN:
        umi_10X = umi + 'C' * (UMI_10X_LEN - len(umi))
    else:
        umi_10X = umi

    seq2_insert = 90
    seq2_cut = 60

    new_seq2_1 = seq2[0:seq2_insert]
    new_seq2_2 = seq2[seq2_cut:] 

    new_seq1 = barcode_10X + umi_10X + TSO
    new_qual1 = 'J' * len(new_seq1)
    new_qual2_1 = qual2[0:len(new_seq2_1)]
    new_qual2_2 = qual2[seq2_cut:]

    return new_seq1, new_qual1, new_seq2_1, new_qual2_1, new_seq2_2, new_qual2_2


def run_convert(BARCODES_10X_FILE, fq2, outdir):

    sample = os.path.abspath(fq2).split('/')[-1]
    sample_prefix = '_'.join(sample.split('_')[:-1])
    barcodes_10X = open(BARCODES_10X_FILE, 'r')
    fq_file = pysam.FastxFile(fq2)

    out_fq1 = xopen(f'{outdir}/{sample_prefix}_S1_L001_R1_001.fastq.gz', 'w')
    out_fq2 = xopen(f'{outdir}/{sample_prefix}_S1_L001_R2_001.fastq.gz', 'w')   
    barcode_dict = {}

    for entry in fq_file:
        name = entry.name
        attrs = name.split('_')
        barcode = attrs[0]
        umi = attrs[1]
        seq = entry.sequence
        qual = entry.quality
        new_seq1, new_qual1, new_seq2_1, new_qual2_1, new_seq2_2, new_qual2_2 = convert_seq(barcode, umi, barcode_dict, barcodes_10X, seq, qual)
        
        out_fq1.write(fastq_line(name, new_seq1, new_qual1))
        out_fq2.write(fastq_line(name, seq, qual))

    out_fq1.close()
    out_fq2.close()
    barcodes_10X.close()

    barcode_record = pd.DataFrame()
    barcode_record['sgr'] = list(barcode_dict.keys())
    barcode_record['10X'] = [barcode_dict[i] for i in barcode_dict]
    barcode_record.to_csv(f'{outdir}/barcode_correspondence_file', sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='convert sgr fq to 10X')
    parser.add_argument('--soft_path', help='cellranger absolute path', required=True)
    parser.add_argument('--sgr_fq', help='read2 fastq file in 01.barcode', required=True)
    parser.add_argument('--outdir', help='output path of 10X format fastq file', required=True)
    args = parser.parse_args()

    check_mkdir(args.outdir)
    BARCODE_10X_FILE = get_10X_barcode(args.soft_path)
    run_convert(BARCODE_10X_FILE, args.sgr_fq, args.outdir)


if __name__ == '__main__':
    main()