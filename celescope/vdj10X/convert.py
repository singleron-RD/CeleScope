from collections import defaultdict

import pandas as pd
import pysam
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from xopen import xopen

BARCODES_10X_FILE = "/SGRNJ/Database/script/soft/cellranger/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/737K-august-2016.txt"
UMI_10X_LEN = 10
TSO = "TTTCTTATATGGG"
SEQ_LEN = 150
BCR_UTR = 'AGTGAGGTCATAGCTCCACCCATTGTAGCTAGCTAGTAGTTTGATTCAGCTCAGCTGTGAGAGAACAGGACCAGGT'


def rev_compl(seq):
    base_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(base_dict[base] for base in reversed(seq))

def convert_seq(sgr_barcode, umi, barcode_dict, barcodes_10X, seq1, seq2, qual2,  method='read2'):
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
    seq1_insert = 14
    if method == 'read2':
        seq1_add = rev_compl(seq2[seq2_insert:])
        new_seq2 = seq2[0:seq2_insert]  
    elif method == 'read1':
        seq1_add = seq1[SEQ_LEN-seq1_insert+1:]
        new_seq2 = seq2
    new_seq1 = barcode_10X + umi_10X + TSO + seq1_add
    new_qual1 = 'J' * len(new_seq1)
    new_qual2 = qual2[0:len(new_seq2)]

    return new_seq1, new_qual1, new_seq2, new_qual2

def fastq_line(name, seq, qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'

class Convert(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # common parameter
        self.outdir = args.outdir
        self.method = args.method
        
        # input
        self.fq2 = args.fq2
        
        # output
        self.out_fq1_file = f'{self.outdir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.out_fq2_file = f'{self.outdir}/{self.sample}_S1_L001_R2_001.fastq.gz'
        self.barcode_correspondence_file = f'{self.outdir}/barcode_correspond.txt'
             
    @utils.add_log
    def run_convert(self):
        # read file
        barcodes_10X = open(BARCODES_10X_FILE, 'r')
        fq_file = pysam.FastxFile(self.fq2)
        
        # open out file
        out_fq1 = xopen(self.out_fq1_file, 'w')
        out_fq2 = xopen(self.out_fq2_file, 'w')
        
        # define var
        barcode_dict = {}

        # write out put fastq file
        for entry in fq_file:
            name = entry.name
            attrs = name.split('_')
            barcode = attrs[0]
            umi = attrs[1]
            seq = entry.sequence
            qual = entry.quality
            seq1 = ''
            new_seq1, new_qual1, new_seq2, new_qual2 = convert_seq(barcode, umi, barcode_dict, barcodes_10X, seq1, seq, qual, method=self.method)
            
            out_fq1.write(fastq_line(name, new_seq1, new_qual1))
            out_fq2.write(fastq_line(name, new_seq2, new_qual2))
            
        out_fq1.close()
        out_fq2.close()
        barcodes_10X.close()
        
        # os.system(f'gzip {out_fq1}')
        # os.system(f'gzip {out_fq2}')
        # write barcode correspond file
        barcode_record = pd.DataFrame()
        barcode_record['sgr'] = list(barcode_dict.keys())
        barcode_record['10X'] = [barcode_dict[i] for i in barcode_dict]
        barcode_record.to_csv(self.barcode_correspondence_file, sep='\t', index=False)
        
    def run(self):
        self.run_convert()
        
def convert(args):
    step_name = 'convert'
    convert_obj = Convert(args, step_name)
    convert_obj.run()
    
def get_opts_convert(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq2', help='R2 read file', required=True)
    parser.add_argument('--method', help='convert method, read1 14bp or read2 60bp', choices=['read1', 'read2'], default='read2')
