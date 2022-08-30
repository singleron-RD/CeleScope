import subprocess
import os

from celescope.tools import utils
from celescope.tools.step import Step, s_common


class Cellranger(Step):
    """
    ## Features
    - Single cell RNA-seq Gene Expression analysis by Cellranger.
    ## Output
    - `03.assemble/{sample}` Cellranger count results.
    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fqs_dir = os.path.abspath(args.fqs_dir)
        self.mem = args.mem
        self.other_param = args.other_param
        self.ref_path = args.ref_path
        self.soft_path = args.soft_path

        # out files
        self.cmd_line = f'{self.outdir}/{self.sample}_cmd_line'

    @utils.add_log
    def cellranger_count(self):
        """Cellranger count"""
        cmd = (
            f'{self.soft_path} count '
            f'--id={self.sample} '
            f'--transcriptome={self.ref_path} '
            f'--fastqs={self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )

        if self.other_param:
            cmd += (' ' + self.other_param)

        self.cellranger_count.logger.info(cmd)
        with open(self.cmd_line, 'w') as f:
            f.write(cmd)
        
        cwd = os.getcwd()
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
        # change dir back to avoid can not find '03.assemble/stat.txt' error
        os.chdir(cwd)

        cmd = (
            f'cp {self.outdir}/{self.sample}/outs/web_summary.html {self.outdir}/..'
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.cellranger_count()


def cellranger(args):
    with Cellranger(args) as runner:
        runner.run()


def get_opts_cellranger(parser, sub_program):
    parser.add_argument('--ref_path', help='reference path for cellranger')
    parser.add_argument('--soft_path', help='soft path for cellranger')
    parser.add_argument('--other_param', help='Other cellranger parameters.', default="")
    parser.add_argument('--mem', help='memory(G)', default=10)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir after convert', required=True)
    return parser