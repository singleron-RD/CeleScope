from celescope.tools import utils
from celescope.tools.step import s_common, Step

from collections import defaultdict
from xopen import xopen
import subprocess
import pysam
import os
import random


BARCODE_10X_LEN = 16
WHITELIST_10X_PATH = [
    "/lib/python/atac/barcodes/737K-cratac-v1.txt.gz",
                      ]


def get_opts_convert(parser, sub_program):
    parser.add_argument('--soft_path', help='soft path for cellranger', required=True)
    parser.add_argument('--gzip', help="Output gzipped fastq files.", action='store_true')
    parser.add_argument('--method', help='bulk or 10X or sgr', required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fq1', help='R2 read file', required=True)
        parser.add_argument('--fq2', help='R2 read file', required=True)
        parser.add_argument('--fq3', help='R2 read file', required=True)
    return parser


class Convert(Step):
    """
    ##Features

    - Convert barcodes to 10X format.

    Output

    - `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

    - `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads as cellranger input.

    - `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads as cellranger input.

    - `02.convert/{sample}_S1_L001_R3_001.fastq.gz` New R3 reads as cellranger input

    - read1, barcode, read2, sample index are associated with R1, R2, R3, I1 respectively.
        R1: Read 1
        R2: Dual index i5 read(10x Barcode)
        R3: Read 2
    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.method = args.method
        self.whitelist_10X_file = os.path.dirname(args.soft_path) + WHITELIST_10X_PATH[0]
        self.whitelist_10X_fh = xopen(self.whitelist_10X_file, 'r')
        self.sgr_tenX = {}

        if self.method == 'bulk':
            self.whitelist_10X_fh = utils.read_one_col(self.whitelist_10X_file)[0][:3000]

        self.out_fq1_file = f'{self.outdir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.out_fq2_file = f'{self.outdir}/{self.sample}_S1_L001_R2_001.fastq.gz'
        self.out_fq3_file = f'{self.outdir}/{self.sample}_S1_L001_R3_001.fastq.gz'
        self.barcode_convert_json = f'{self.outdir}/barcode_convert.json'
    
    @utils.add_log
    def write_fq1(self):
        cmd = f'cp {self.args.fq1} {self.out_fq1_file}'
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def write_fq2(self):
        if self.method == '10X':
            cmd = f'cp {self.args.fq2} {self.out_fq2_file}'
            subprocess.check_call(cmd, shell=True)

        else:
            out_fq2 = xopen(self.out_fq2_file, 'w')
            with pysam.FastxFile(self.args.fq2) as fq2_fh:
                for entry in fq2_fh:
                    name = entry.name
                    comment = entry.comment
                    sgr_barcode = entry.sequence

                    if self.method=='bulk':
                        barcode_10X = random.choice(self.whitelist_10X_fh)
                    else:
                        if sgr_barcode in self.sgr_tenX:
                            barcode_10X = self.sgr_tenX[sgr_barcode]
                        else:
                            # new barcode from whitelist
                            barcode_10X = self.whitelist_10X_fh.readline().strip()
                            self.sgr_tenX[sgr_barcode] = barcode_10X

                    new_name = f'{name} {comment}'
                    new_qual = 'F' * len(barcode_10X)
                    out_fq2.write(f'@{new_name}\n{barcode_10X}\n+\n{new_qual}\n')

            out_fq2.close()
    
    @utils.add_log
    def write_fq3(self):
        cmd = f'cp {self.args.fq3} {self.out_fq3_file}'
        subprocess.check_call(cmd, shell=True)
    
    @utils.add_log
    def dump_tenX_sgr_barcode_json(self):
        tenX_sgr = {}
        for sgr, tenX in self.sgr_tenX.items():
            if self.bulk_seq:
                tenX = ','.join(tenX)
            tenX_sgr[tenX] = sgr

        utils.dump_dict_to_json(tenX_sgr, self.barcode_convert_json)
    
    def run(self):
        self.write_fq1()
        self.write_fq2()
        self.write_fq3()
        if self.method == 'sgr':
            self.dump_tenX_sgr_barcode_json()

def convert(args):
    step_name = 'convert'
    with Convert(args, step_name) as runner:
        runner.run()

