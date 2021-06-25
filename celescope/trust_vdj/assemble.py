import os
from celescope.tools import utils
from celescope.tools.Step import Step, s_common
import pandas as pd


TRUST = '/SGRNJ03/randd/zhouxin/software/TRUST4/'


class Assemble(Step):
    """
    Features

    - Assemble TCR/BCR
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.sample = args.sample
        self.species = args.species
        self.speed_up = args.speed_up


    @utils.add_log
    def run(self):

        species = self.species

        index_file = f'{TRUST}/index/{species}/{species}_ref.fa'
        ref = f'{TRUST}/index/{species}/{species}_IMGT+C.fa'

        string1 = ''
        if self.speed_up:
            string1 = '--repseq '
        cmd = (
            f'{TRUST}/run-trust4 -t {self.thread} '
            f'-u {self.fq2} '
            f'--barcode {self.fq1} '
            f'--barcodeRange 0 23 + '
            f'-f {index_file} '
            f'--ref {ref} '
            f'{string1}'
            f'-o {self.sample} --od {self.outdir}/TRUST4' 
        )

        Assemble.run.logger.info(cmd)

        if not os.path.exists(f'{self.outdir}/TRUST4/{self.sample}_barcode_report.tsv'):
            os.system(cmd)

            #fq = f'{self.outdir}/TRUST4/{self.sample}_toassemble.fq'


@utils.add_log
def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='R1 reads from match step', required=True)
        parser.add_argument('--fq2', help='R2 reads from match step', required=True)

    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)
    parser.add_argument('--speed_up', help='speed assemble for TCR/BCR seq data', action='store_true')       








