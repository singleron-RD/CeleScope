"""
split scRNA-Seq fastq file(01.barcode/{sample}_2.fq)
"""
import glob
import os
from collections import defaultdict

import pysam
import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT
from celescope.tools.cellranger3.wrapper import Cell_calling, read_raw_matrix


class Split_tag(Step):
    """
    Features
    - Split scRNA-Seq fastq according to tag assignment.

    Output
    - `matrix/` Matrix files of each tag.(Optional)
    - `fastq/` Fastq files of each tag.(Optional)
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        if not (args.split_matrix or args.split_fastq):
            return

        # set
        df_umi_tag = pd.read_csv(args.umi_tag_file, sep='\t', index_col=0)
        df_umi_tag = df_umi_tag.rename_axis('barcode').reset_index()
        self.tag_barcode_dict = {tag: set(row["barcode"].tolist()) for tag, row in df_umi_tag.groupby("tag")}

        if args.split_matrix:
            self.matrix_outdir = f'{args.outdir}/matrix/'
            if args.match_dir:
                matrix_10X_dir = glob.glob(f'{args.match_dir}/05.count/*_matrix_10X*')[0]
            elif args.matrix_dir:
                matrix_10X_dir = args.matrix_dir
            else:
                raise ValueError("--match_dir or --matrix_dir is required.")
            self.raw_mat, self.raw_features_path, self.raw_barcodes = read_raw_matrix(matrix_10X_dir)

        if args.split_fastq:
            self.rna_fq_file = glob.glob(f'{args.match_dir}/*barcode/*_2.fq*')[0]

            fastq_outdir = f'{args.outdir}/fastqs/'
            os.system(f'mkdir -p {fastq_outdir}')

            self.r2_fastq_files_handle = {}
            self.r1_fastq_files_handle = {}
            for tag in self.tag_barcode_dict:
                r2_fastq_file_name = f'{fastq_outdir}/{tag}_2.fq'
                self.r2_fastq_files_handle[tag] = open(r2_fastq_file_name, 'w')
                r1_fastq_file_name = f'{fastq_outdir}/{tag}_1.fq'
                self.r1_fastq_files_handle[tag] = open(r1_fastq_file_name, 'w')

            self.tag_read_index_dict = defaultdict(set)

    @utils.add_log
    def write_r2_fastq_files(self):
        read_num = 0
        with pysam.FastxFile(self.rna_fq_file, 'r') as rna_fq:
            for read in rna_fq:
                read_num += 1
                attr = read.name.strip("@").split("_")
                barcode = attr[0]
                read_index = int(attr[2])
                for tag in self.tag_barcode_dict:
                    if barcode in self.tag_barcode_dict[tag]:
                        self.tag_read_index_dict[tag].add(read_index)
                        self.r2_fastq_files_handle[tag].write(str(read) + '\n')

                if read_num % 1000000 == 0:
                    self.write_r2_fastq_files.logger.info(f'{read_num} done')

        for tag in self.r2_fastq_files_handle:
            self.r2_fastq_files_handle[tag].close()

    @utils.add_log
    def write_r1_fastq_files(self):
        with pysam.FastxFile(self.args.R1_read, 'r') as r1_read:
            for read_index, read in enumerate(r1_read, start=1):
                for tag in self.tag_read_index_dict:
                    if read_index in self.tag_read_index_dict[tag]:
                        self.r1_fastq_files_handle[tag].write(str(read) + '\n')

        for tag in self.r1_fastq_files_handle:
            self.r1_fastq_files_handle[tag].close()

    @utils.add_log
    def split_matrix(self):
        for tag in self.tag_barcode_dict:
            outdir = f'{self.matrix_outdir}/{tag}_matrix_10X/'
            runner = Cell_calling(outdir, self.raw_mat, self.raw_features_path, self.raw_barcodes)
            tag_barcodes = list(self.tag_barcode_dict[tag])
            raw_barcodes = list(runner.raw_barcodes)
            tag_barcodes_indices = [raw_barcodes.index(barcode) for barcode in tag_barcodes]
            tag_barcodes_indices.sort()
            runner.write_slice_matrix(tag_barcodes_indices)


    @utils.add_log
    def run(self):
        if self.args.split_matrix:
            self.split_matrix()
        if self.args.split_fastq:
            self.write_r2_fastq_files()
            self.write_r1_fastq_files()


def split_tag(args):
    step_name = "split_tag"
    runner = Split_tag(args, step_name)
    runner.run()


def get_opts_split_tag(parser, sub_program):
    parser.add_argument(
        "--split_fastq",
        help="If used, will split scRNA-Seq fastq file according to tag assignment.",
        action='store_true',
    )
    parser.add_argument(
        "--split_matrix",
        help="If used, will split scRNA-Seq matrix file according to tag assignment.",
        action='store_true',
    )
    if sub_program:
        parser.add_argument("--umi_tag_file", help="UMI tag file.", required=True)
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--matrix_dir", help="Match celescope scRNA-Seq matrix directory.")
        parser.add_argument("--R1_read", help='R1 read path.')
        s_common(parser)
