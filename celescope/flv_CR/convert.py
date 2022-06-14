import os
import pysam
import pandas as pd
from xopen import xopen

from celescope.tools import utils
from celescope.tools.step import Step, s_common


UMI_10X_LEN = 10
TSO = "TTTCTTATATGGG"


class Convert(Step):
    """
    Features

    - Convert barcodes and UMI to 10X format.

    Output        
    - `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

    - `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads as cellranger input.

    - `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads as cellranger input.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.fq2 = args.fq2
        self.not_split_R2 = args.not_split_R2

        if args.soft_path:
            self.soft_path = args.soft_path
        else:
            self.soft_path = f'/SGRNJ/Database/script/soft/cellranger/cellranger-{args.version}/cellranger'
        self.whitelist_10X = os.path.dirname(self.soft_path) + "/lib/python/cellranger/barcodes/737K-august-2016.txt"

        # out
        self.out_fq1_file = f'{self.outdir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.out_fq2_file = f'{self.outdir}/{self.sample}_S1_L001_R2_001.fastq.gz'

    @utils.add_log
    def run_convert(self):
        """Convert fastq file to new file in 10X format"""
        
        barcode_dict = {}
        whitelist_10X = open(self.whitelist_10X, 'r')
        fq_file = pysam.FastxFile(self.fq2)
        
        out_fq1 = xopen(self.out_fq1_file, 'w')
        out_fq2 = xopen(self.out_fq2_file, 'w')
        
        for entry in fq_file:
            name = entry.name
            attrs = name.split('_')
            barcode, umi = attrs[0], attrs[1]
            seq, qual = entry.sequence, entry.quality
            new_seq1, new_qual1, new_seq2_1, new_qual2_1, new_seq2_2, new_qual2_2 = Convert.convert_seq(barcode, umi, barcode_dict, whitelist_10X, seq, qual)

            if not self.not_split_R2:
                out_fq1.write(f'@{name}_1\n{new_seq1}\n+\n{new_qual1}\n')
                out_fq1.write(f'@{name}_2\n{new_seq1}\n+\n{new_qual1}\n')
                out_fq2.write(f'@{name}_1\n{new_seq2_1}\n+\n{new_qual2_1}\n')
                out_fq2.write(f'@{name}_2\n{new_seq2_2}\n+\n{new_qual2_2}\n')
            else:
                out_fq1.write(f'@{name}\n{new_seq1}\n+\n{new_qual1}\n')
                out_fq2.write(f'@{name}_1\n{seq}\n+\n{qual}\n')

        out_fq1.close()
        out_fq2.close()
        whitelist_10X.close()

        barcode_record = pd.DataFrame()
        barcode_record['sgr'] = list(barcode_dict.keys())
        barcode_record['10X'] = [barcode_dict[i] for i in barcode_dict]
        barcode_record.to_csv(f'{self.outdir}/barcode_correspond.txt', sep='\t', index=False)
    
    @staticmethod
    def convert_seq(sgr_barcode, sgr_umi, barcode_dict, barcodes_10X, seq2, qual2):
        """Convert each seq to 10X format

        :param sgr_barcode: sgr barcode
        :param sgr_umi: sgr umi
        :param barcode_dict: key:SGR barcode; value:10X barcode
        :param barcodes_10X: 10X barcode
        :param seq2: r2 sequence
        :param qual2: r2 sequence quality
        :return: sequence in 10X format
        """

        if sgr_barcode in barcode_dict:
            barcode_10X = barcode_dict[sgr_barcode]
        else:
            barcode_10X = barcodes_10X.readline().strip()
            barcode_dict[sgr_barcode] = barcode_10X

        if len(sgr_umi) > UMI_10X_LEN:
            umi_10X = sgr_umi[:UMI_10X_LEN]
        elif len(sgr_umi) < UMI_10X_LEN:
            umi_10X = sgr_umi + 'C' * (UMI_10X_LEN - len(sgr_umi))
        else:
            umi_10X = sgr_umi

        new_seq2_1, new_seq2_2 = seq2[:90], seq2[60:]

        new_seq1 = barcode_10X + umi_10X + TSO
        new_qual1 = 'J' * len(new_seq1)
        new_qual2_1 = qual2[:len(new_seq2_1)]
        new_qual2_2 = qual2[60:]

        return new_seq1, new_qual1, new_seq2_1, new_qual2_1, new_seq2_2, new_qual2_2

    def run(self):
        self.run_convert()


def convert(args):
    step_name = 'convert'
    convert_obj = Convert(args, step_name)
    convert_obj.run()


def get_opts_convert(parser, sub_program):
    parser.add_argument('--soft_path', help='soft path for cellranger')
    parser.add_argument('--not_split_R2', help='whether to split R2',action='store_true')
    parser.add_argument('--version', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'],
                        default='4.0.0')
    if sub_program:
        s_common(parser)
        parser.add_argument('--fq2', help='R2 read file', required=True)
    return parser
