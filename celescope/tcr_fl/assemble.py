
import os
from os import listdir
from os.path import isfile, join
from concurrent.futures import ProcessPoolExecutor

from celescope.tools.step import Step, s_common
from celescope.tools.utils import add_log

TRACER_PATH = '/SGRNJ/Public/Software/tracer/tracer'
CONF_PATH = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/tcr_fl/20201103/tracer_SGR.conf'
CONDA = 'vdjpuzzle1'
CONDA_SUB = 'celescope_tracer'

@add_log
def tracer_summarise(outdir):
    tracer_outdir = f'{outdir}/tracer'
    cmd = (
        f'source activate {CONDA_SUB}; '
        f'{TRACER_PATH} summarise '
        f'-c {CONF_PATH} '
        f'{tracer_outdir} '
    )
    tracer_summarise.logger.info(cmd)
    os.system(cmd)


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
    os.system(cmd)
    

@add_log
def run_assemble(outdir, fastq_dir, thread):
    fqs = [join(fastq_dir, f) for f in listdir(fastq_dir) if isfile(join(fastq_dir, f))]
    outdirs = [outdir] * len(fqs)
    if not os.path.exists(f'{outdir}/tracer'):
        os.makedirs(f'{outdir}/tracer')
    
    all_res = []
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(tracer, fqs, outdirs):
            all_res.append(res)
    
    tracer_summarise(outdir)


def assemble(args):
    thread = int(args.thread)
    fastq_dir = args.fastq_dir
    outdir = args.outdir
    run_assemble(outdir, fastq_dir, thread)


def get_opts_assemble(parser, sub_program):
    s_common(parser)
    if sub_program:
        parser.add_argument("--fastq_dir", required=True)

