"""
split scRNA-Seq fastq file(01.barcode/{sample}_2.fq)
"""
import glob
import os

import pysam
import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT

class Split_tag(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # set

        df_umi_tag = pd.read_csv(args.umi_tag_file, sep='\t', index_col=0)
        df_umi_tag = df_umi_tag.rename_axis('barcode').reset_index()
        self.tag_barcode_dict = {tag: set(row["barcode"].tolist()) for tag, row in df_umi_tag.groupby("tag")}

        if args.split_fastq:
            self.rna_fq_file = glob.glob(f'{args.match_dir}/*barcode/*_2.fq*')[0]

            fastq_outdir = f'{args.outdir}/fastqs/'
            os.system(f'mkdir -p {fastq_outdir}')
            self.fastq_files_handle = {}
            for tag in self.tag_barcode_dict:
                fastq_file_name = f'{fastq_outdir}/{tag}_2.fq'
                self.fastq_files_handle[tag] = open(fastq_file_name, 'w')

    @utils.add_log
    def write_fastq_files(self):
        read_num = 0
        with pysam.FastxFile(self.rna_fq_file, 'r') as rna_fq:
            for read in rna_fq:
                read_num += 1
                attr = read.name.strip("@").split("_")
                barcode = attr[0]
                for tag in self.tag_barcode_dict:
                    if barcode in self.tag_barcode_dict[tag]:
                        self.fastq_files_handle[tag].write(str(read))

                if read_num % 1000000 == 0:
                    self.write_fastq_files.logger.info(f'{read_num} done')

        for tag in self.tag_barcode_dict:
            self.fastq_files_handle[tag].close()


    @utils.add_log
    def run(self):
        if self.args.split_fastq:
            self.write_fastq_files()

def split_tag(args):
    step_name = "split_tag"
    runner = Split_tag(args, step_name)
    runner.run()

def get_opts_split_tag(parser, sub_program):
    parser.add_argument(
        "--split_fastq", 
        help="Split scRNA-Seq fastq file(01.barcode/{sample}_2.fq).",
        action='store_true',
    )
    if sub_program:
        parser.add_argument("--umi_tag_file", help="UMI tag file", required=True)
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'], required=True)
        s_common(parser)

