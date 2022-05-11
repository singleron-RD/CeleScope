import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from os import listdir
from os.path import isfile, join

from celescope.tools.step import Step, s_common
from celescope.tools.utils import add_log

TRACER_PATH = '/SGRNJ/Public/Software/tracer/tracer'
CONF_PATH = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/tcr_fl/20201103/tracer_SGR.conf'
CONDA = 'vdjpuzzle1'
CONDA_SUB = 'celescope_tracer'


def tracer(fq, outdir):
    prefix = os.path.basename(fq).strip('.fq')
    cmd = (
        f'source activate {CONDA}; '
        f'{TRACER_PATH} assemble '
        f'--fragment_length 150 '
        f'--fragment_sd 5 '
        f'--single_end '
        f'--species Hsap '
        f'-c {CONF_PATH} '
        f'{fq} '
        f'{prefix} '
        f'{outdir}/tracer '
    )
    subprocess.check_call(cmd, shell=True)


class Assemble_TCR(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.fastq_dir = args.fastq_dir

    @add_log
    def tracer_summarise(self):
        tracer_outdir = f'{self.outdir}/tracer'
        cmd = (
            f'source activate {CONDA_SUB}; '
            f'{TRACER_PATH} summarise '
            f'-c {CONF_PATH} '
            f'{tracer_outdir} '
        )
        Assemble_TCR.tracer_summarise.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @add_log
    def run(self):
        fqs = [join(self.fastq_dir, f) for f in listdir(self.fastq_dir) if isfile(join(self.fastq_dir, f))]
        outdirs = [self.outdir] * len(fqs)
        if not os.path.exists(f'{self.outdir}/tracer'):
            os.makedirs(f'{self.outdir}/tracer')

        all_res = []
        with ProcessPoolExecutor(self.thread) as pool:
            for res in pool.map(tracer, fqs, outdirs):
                all_res.append(res)
        self.tracer_summarise()
        self._clean_up()


def assemble(args):
    with Assemble_TCR(args) as runner:
        runner.run()


def get_opts_assemble(parser, sub_program):
    s_common(parser)
    if sub_program:
        parser.add_argument("--fastq_dir", required=True)
