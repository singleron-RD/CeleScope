import os
import subprocess
import glob
import pandas as pd

from celescope.tools.step import Step, s_common
from celescope.tools import utils


class Assemble(Step):
    """
    ## Features

    - TCR/BCR Assemble by Cellranger.

    ## Output
    
    - `03.assemble/{sample}` Cellranger vdj results.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fqs_dir = os.path.abspath(args.fqs_dir)
        self.mem = args.mem
        self.other_param = args.other_param
        if args.ref_path and args.soft_path:
            self.ref_path = args.ref_path
            self.soft_path = args.soft_path
        else:
            self.soft_path = f'/SGRNJ/Database/script/soft/cellranger/cellranger-{args.version}/cellranger'
            self.ref_path = glob.glob(f'/SGRNJ/Database/script/soft/cellranger/vdj_ref/{args.version}/{args.species}/refdata*')
            if not self.ref_path:
                self.ref_path = glob.glob(f'/SGRNJ/Database/script/soft/cellranger/vdj_ref/4.0.0/{args.species}/refdata*')
            self.ref_path = [i for i in self.ref_path if os.path.isdir(i)][0]

        self.not_split_R2 = args.not_split_R2
        self.seqtype = args.seqtype
        if self.seqtype == 'TCR':
            self.chains = ['TRA', 'TRB']
            self.pair = ['TRA_TRB']
        elif self.seqtype == 'BCR':
            self.chains = ['IGH', 'IGL', 'IGK']

    @utils.add_log
    def run_assemble(self):
        """Cellranger vdj"""
        cmd = (
            f'{self.soft_path} vdj '
            f'--id={self.sample} '
            f'--reference={self.ref_path} '
            f'--fastqs={self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )

        if self.other_param:
            cmd += (self.other_param)

        self.run_assemble.logger.info(cmd)
        with open(f'{self.outdir}/{self.sample}_cmd_line', 'w') as f:
            f.write(cmd)
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
    
    @utils.add_log
    def gen_report(self):
        stat_dict = pd.read_csv(f'{self.outdir}/../01.barcode/stat.txt', sep=':', header=None)
        read_count = int(stat_dict.iloc[0, 1].replace(',', ''))
        sum_dict = pd.read_csv(f'{self.outdir}/{self.sample}/outs/metrics_summary.csv', sep=',', index_col=None)
        sum_dict = sum_dict.T.to_dict()
        total_reads = int(sum_dict[0]["Number of Read Pairs"].replace(',', ''))
        cell_nums = len(set(self.filter_contig.barcode))

        _index = 200
        if self.not_split_R2:
            _index = 100

        self.add_metric(
            name='Estimated Number of Cells',
            value=cell_nums,
            help_info=f"Cells with at least one productive {' or '.join(self.seqtype)} chain"
        )

        self.add_metric(
            name='Reads Mapped to Any V(D)J Gene',
            value=int(total_reads * (float(sum_dict[0]['Reads Mapped to Any V(D)J Gene'].strip('%'))/_index)),
            total=read_count,
            help_info=f"Reads mapped to any {' or '.join(self.seqtype)} genes."
        )

        for chain in self.chains:
            self.add_metric(
                name=f'Reads Mapped to {chain}',
                value=int(
                    total_reads * (float(sum_dict[0][f'Reads Mapped to {chain}'].strip('%'))/_index)),
                total=read_count,
                help_info=f"Reads mapped to {chain} chain. For BCR, this should be one of {self.seqtype}"
            )

        self.add_metric(
            name='Fraction Reads in Cells',
            value=int(total_reads * (float(sum_dict[0]['Fraction Reads in Cells'].strip('%'))/_index)),
            total=read_count,
            help_info="Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes"
        )

        for chain in self.chains:
            mid = self.filter_contig[self.filter_contig['chain']== chain]['umis'].median()
            if mid == mid:
                self.add_metric(
                    name=f'Median used {chain} UMIs per Cell',
                    value=int(mid),
                    help_info=f"Median number of UMIs assigned to a {chain} contig per cell."
                )
            else:
                self.add_metric(
                    name=f'Median used {chain} UMIs per Cell',
                    value=0,
                    help_info=f"Median number of UMIs assigned to a {chain} contig per cell."
                )

    def run(self):
        self.run_assemble()


def assemble(args):
    assemble_obj = Assemble(args)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    parser.add_argument('--species', help='species', choices=['hs', 'mmu'], default='hs', required=True)
    parser.add_argument('--version', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'],
                        default='4.0.0')
    parser.add_argument('--ref_path', help='reference path for cellranger')
    parser.add_argument('--soft_path', help='soft path for cellranger')
    parser.add_argument('--other_param', help='Other cellranger parameters.', default="")
    parser.add_argument('--mem', help='memory(G)', default=10)
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--not_split_R2', help='whether split r2',action='store_true')
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir after convert', required=True)
    return parser

