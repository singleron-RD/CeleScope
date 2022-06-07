import glob
import gzip
import argparse
import os

import pysam

from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.celescope import ArgFormatter


SAMPLE = 'test1'

class Extract_read():
    """
    Extract fq1 and fq2 reads according to barcode.
    """

    def __init__(self, args):
        self.args = args
        self.rna_fq_file = glob.glob(f'{args.match_dir}/*barcode/*_2.fq*')[0]
        self.barcodes, _num = utils.read_one_col(args.barcode_file)
        self.read_index = set()

        # out
        self.out_fq1 = f'{args.outdir}/{SAMPLE}_1.fq.gz'
        self.out_fq2 = f'{args.outdir}/{SAMPLE}_2.fq.gz'

        # mkdir
        if not os.path.exists(args.outdir):
            os.system(f'mkdir -p {args.outdir}')


    @utils.add_log
    def write_r2_fastq_files(self):
        read_num = 0
        with pysam.FastxFile(self.rna_fq_file, 'r') as rna_fq:
            with gzip.open(self.out_fq2, 'wt') as fq2:
                for read in rna_fq:
                    read_num += 1
                    attr = read.name.strip("@").split("_")
                    barcode = attr[0]
                    read_index = int(attr[2])
                    if barcode in self.barcodes:
                        self.read_index.add(read_index)
                        fq2.write(str(read) + '\n')

                    if read_num % 1000000 == 0:
                        self.write_r2_fastq_files.logger.info(f'{read_num} done')


    @utils.add_log
    def write_r1_fastq_files(self):
        with pysam.FastxFile(self.args.R1_read, 'r') as r1_read:
            with gzip.open(self.out_fq1, 'wt') as fq1:
                for read_index, read in enumerate(r1_read, start=1):
                    if read_index in self.read_index:
                        fq1.write(str(read) + '\n')

    def run(self):
        self.write_r2_fastq_files()
        self.write_r1_fastq_files()


def main():
    parser = argparse.ArgumentParser(description='Extract fq1 and fq2 reads according to barcode.', formatter_class=ArgFormatter) 
    parser.add_argument("--barcode_file", help="file with barcode to extract", required=True)
    parser.add_argument("--match_dir", help=HELP_DICT['match_dir'], required=True)
    parser.add_argument("--R1_read", help='R1 read path.')
    parser.add_argument("--outdir", help=HELP_DICT['outdir'], default='./')

    args = parser.parse_args()
    runner = Extract_read(args)
    runner.run()

if __name__ == '__main__':
    main()