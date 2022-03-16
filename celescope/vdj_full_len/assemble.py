import os
import subprocess

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj_full_len.__init__ import ref_dict, soft_dict


class Assemble(Step):
    """
    ## Features

    - TCR/BCR Assemble.

    ## Output
    - `03.assemble/{sample}/outs/` Recording assemble results.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # common parameters
        self.mem = args.mem
        self.species = args.species
        self.soft = args.soft
        self.thread = args.thread
        self.fqs_dir = args.fqs_dir
        self.other_param = args.other_param

        self.cwd_path = os.getcwd()

        self.ref_path = ref_dict[self.soft][self.species]
        self.soft_path = soft_dict[self.soft]
        
        # for customers
        if args.ref_path and args.soft_path:
            self.ref_path = args.ref_path
            self.soft_path = args.soft_path

    @utils.add_log
    def run_assemble(self):

        cmd = (
            f'{self.soft_path} vdj '
            f'--id={self.sample} '
            f'--reference={self.ref_path} '
            f'--fastqs={self.cwd_path}/{self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )

        if self.other_param:
            cmd += (self.other_param)

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
                        choices=['hs', 'mmu'], default='hs')
    parser.add_argument('--soft', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'],
                        default='4.0.0')
    parser.add_argument('--mem', help='memory (G)', default=10)
    parser.add_argument('--ref_path', help='reference path for cellranger')
    parser.add_argument('--soft_path', help='soft path for cellranger')
    parser.add_argument('--other_param', help='Other cellranger parameters.', default="")
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir', required=True)
    return parser
