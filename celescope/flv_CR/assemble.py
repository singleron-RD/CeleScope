import os
import subprocess
import glob
import pandas as pd

from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.flv_trust4.__init__ import CHAIN, PAIRED_CHAIN


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
        self.ref_path = args.ref_path
        self.soft_path = args.soft_path
        self.seqtype = args.seqtype

        self.chains = CHAIN[self.seqtype]
        self.pair = PAIRED_CHAIN[self.seqtype]

        # out files
        self.cmd_line = f'{self.outdir}/{self.sample}_cmd_line'


    @utils.add_log
    def assemble(self):
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

        self.assemble.logger.info(cmd)
        with open(self.cmd_line, 'w') as f:
            f.write(cmd)
        
        cwd = os.getcwd()
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
        # change dir back to avoid can not find '03.assemble/stat.txt' error
        os.chdir(cwd)

    def run(self):
        self.assemble()

def assemble(args):
    with Assemble(args) as runner:
        runner.run()

    with Mapping(args) as runner:
        runner.run()

    # TODO
    with Cells(args) as runner:
        runner.run()

    with Annotation(args) as runner:
        runner.run()


class cellranger_metrics(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        cellranger_metrics_csv = f'{self.outdir}/{self.sample}/outs/metrics_summary.csv'
        df = pd.read_csv(cellranger_metrics_csv, index_col=None)
        self.metrics_dict = df.T.to_dict()[0]

        self.seqtype = args.seqtype
        self.chains = CHAIN[self.seqtype]

    def run(self):
        self.add_metrics()


class Mapping(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)


    def add_metrics(self):
        name = 'Reads Mapped to Any V(D)J Gene'
        self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info="Fraction of reads that partially or wholly map to any germline V(D)J gene segment."
        )

        for chain in self.chains:
            name = f'Reads Mapped to {chain}'
            self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info=f"Fraction of reads that map partially or wholly to a germline {chain} gene segment."
        )

class Cells(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    def add_metrics(self):
        pass

class Annotation(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
    
    def add_metrics(self):
        pass


    
def get_opts_assemble(parser, sub_program):
    parser.add_argument('--ref_path', help='reference path for cellranger')
    parser.add_argument('--soft_path', help='soft path for cellranger')
    parser.add_argument('--other_param', help='Other cellranger parameters.', default="")
    parser.add_argument('--mem', help='memory(G)', default=10)
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir after convert', required=True)
    return parser

