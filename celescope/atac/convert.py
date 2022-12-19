from celescope.tools import utils
from celescope.tools.step import s_common, Step

from xopen import xopen
import subprocess
import pysam
import os


# Const
BARCODE_10X_LEN = 16
WHITELIST_10X_PATH = "/lib/python/atac/barcodes/737K-cratac-v1.txt.gz"


def get_opts_convert(parser, sub_program):
    parser.add_argument('--soft_path', help='soft path for cellranger', required=True)
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
    - `02.convert/barcode_correspond.txt` Recording barcodes correspondence.(if --method sgr)
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

        self.whitelist_10X_file = os.path.dirname(args.soft_path) + WHITELIST_10X_PATH

        self.whitelist_10X_fh = xopen(self.whitelist_10X_file, 'r')
        self.sgr_tenX = {}

        # self.whitelist_10X_fh = iter(utils.read_one_col(self.whitelist_10X_file))

        self.out_fq1_file = f'{self.outdir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.out_fq2_file = f'{self.outdir}/{self.sample}_S1_L001_R2_001.fastq.gz'
        self.out_fq3_file = f'{self.outdir}/{self.sample}_S1_L001_R3_001.fastq.gz'
        self.barcode_convert_json = f'{self.outdir}/barcode_convert.json'
    
    @utils.add_log
    def gzip_fq1(self):
        cmd = f'gzip -c {self.args.fq1} > {self.out_fq1_file}'
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def write_fq2(self):
        out_fq2 = xopen(self.out_fq2_file, 'w')

        with pysam.FastxFile(self.args.fq2) as fq2_fh:
            for entry in fq2_fh:
                name = entry.name
                sgr_barcode = entry.sequence

                if sgr_barcode in self.sgr_tenX:
                    barcode_10X = self.sgr_tenX[sgr_barcode]
                else:
                    # new barcode from whitelist
                    barcode_10X = self.whitelist_10X_fh.readline().strip()
                    self.sgr_tenX[sgr_barcode] = barcode_10X

                new_qual = 'F' * len(barcode_10X)
                out_fq2.write(f'@{name}\n{barcode_10X}\n+\n{new_qual}\n')

            out_fq2.close()
    
    @utils.add_log
    def gzip_fq3(self):
        cmd = f'gzip -c {self.args.fq3} > {self.out_fq3_file}'
        subprocess.check_call(cmd, shell=True)
    
    @utils.add_log
    def dump_tenX_sgr_barcode_json(self):
        tenX_sgr = {}
        for sgr, tenX in self.sgr_tenX.items():
            tenX_sgr[tenX] = sgr

        utils.dump_dict_to_json(tenX_sgr, self.barcode_convert_json)
    
    def run(self):
        self.gzip_fq1()
        self.write_fq2()
        self.gzip_fq3()
        self.dump_tenX_sgr_barcode_json()


def convert(args):
    step_name = 'convert'
    with Convert(args, step_name) as runner:
        runner.run()