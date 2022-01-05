import os
import subprocess

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj_full_len.__init__ import ref_dict, soft_dict


class Assemble(Step):
    """
    run CellRanger 
    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # common parameters
        self.mem = args.mem
        self.species = args.species
        self.soft = args.soft
        self.thread = args.thread
        self.fqs_dir = args.fqs_dir

        self.cwd_path = os.getcwd()

    @utils.add_log
    def run_assemble(self):
        ref_path = ref_dict[self.soft][self.species]
        soft_path = soft_dict[self.soft]
        cmd = (
            f'{soft_path} vdj '
            f'--id={self.sample} '
            f'--reference={ref_path} '
            f'--fastqs={self.cwd_path}/{self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )
        Assemble.run_assemble.logger.info(cmd)
        with open(f'{self.outdir}/{self.sample}_vdj_10X.sh', 'w') as f:
            f.write(cmd)
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_assemble()

def assemble(args):
    assemble_obj = Assemble(args)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    parser.add_argument('--species', help='species',
                        choices=['hs', 'mmu'], required=True)
    parser.add_argument('--soft', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'],
                        default='4.0.0')
    parser.add_argument('--mem', help='memory (G)', default=10)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir', required=True)
    return parser
