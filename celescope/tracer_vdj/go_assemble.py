import argparse
import os
from os import listdir
from os.path import isfile, join
from concurrent.futures import ProcessPoolExecutor
from celescope.tools.utils import add_log
import datetime


TRACER_PATH = '/SGRNJ03/randd/zhouxin/software/tracer/tracer'
CONF_PATH = '/SGRNJ03/randd/zhouxin/software/tracer/tracer.conf'
BRACER_PATH = '/SGRNJ03/randd/zhouxin/software/bracer/bracer'
BRACER_CONDA = 'bracer'
BRACER_CONF = '/SGRNJ03/randd/zhouxin/software/bracer/bracer.conf'


# 开始组装


def bracer_summarise(outdir):
    bracer_outdir = f'{outdir}/bracer'
    cmd = (
        f'source activate {BRACER_CONDA}; '
        f'{BRACER_PATH} summarise '
        f'-c {BRACER_CONF} '
        f'--no_networks ' 
        f'{bracer_outdir} '
        )
    bracer_summarise.logger.info(cmd)
    os.system(cmd)


def bracer(fq, outdir, species):
    prefix = os.path.basename(fq).strip('.fq')
    cmd = (
        f'source activate {BRACER_CONDA}; '
        f'{BRACER_PATH} assemble '
        f'--fragment_length 150 '
        f'--fragment_sd 5 '
        f'--single_end '
        f'--small_index '
        f'--species {species} '
        f'-c {BRACER_CONF} '
        f'{prefix} '
        f'{outdir}/bracer '
        f'{fq} '
    )
    bracer.logger.info(cmd)
    os.system(cmd)


def tracer_summarise(outdir):
    tracer_outdir = f'{outdir}/tracer'
    cmd = (
        f'source activate {BRACER_CONDA}; '
        f'{TRACER_PATH} summarise '
        f'-c {CONF_PATH} '
        f'--no_networks '
        f'{tracer_outdir} '
    )
    tracer_summarise.logger.info(cmd)
    os.system(cmd)


def tracer(fq, outdir, species):
    prefix = os.path.basename(fq).strip('.fq')
    cmd = (
        f'source activate {BRACER_CONDA}; '
        f'{TRACER_PATH} assemble '
        f'--fragment_length 150 '
        f'--fragment_sd 5 '
        f'--single_end '
        f'--small_index '
        f'-m assembly '
        f'--species {species} '
        f'-c {CONF_PATH} '
        f'{fq} '
        f'{prefix} '
        f'{outdir}/tracer '
    )
    tracer.logger.info(cmd)
    os.system(cmd)


@add_log
def run_tracer(outdir, fastq_dir, species, thread):

    fqs = [join(fastq_dir, f) for f in listdir(fastq_dir) if isfile(join(fastq_dir, f))]
    outdirs = [outdir] * len(fqs)
    species = [species] * len(fqs)
    if not os.path.exists(f'{outdir}/tracer'):
        os.makedirs(f'{outdir}/tracer')

    all_res = []
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(tracer, fqs, outdirs, species):
            all_res.append(res)

    tracer_summarise(outdir)


@add_log
def run_bracer(outdir, fastq_dir, species, thread):
    fqs = [join(fastq_dir, f) for f in listdir(fastq_dir) if isfile(join(fastq_dir, f))]
    outdirs = [outdir] * len(fqs)
    species = [species] * len(fqs)
    if not os.path.exists(f'{outdir}/bracer'):
        os.makedirs(f'{outdir}/bracer')

    all_res = []
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(bracer, fqs, outdirs, species):
            all_res.append(res)

    bracer_summarise(outdir)


def go_assemble(args):
    thread = int(args.thread)
    fastq_dir = args.fastq_dir
    outdir = args.outdir
    species = args.species
    
    mode = args.mode
    if mode == 'TCR':
        run_tracer(outdir, fastq_dir, species, thread)
    elif mode == 'BCR':
        run_bracer(outdir, fastq_dir, species, thread)


def get_opts_go_assemble(parser, sub_program):
    if sub_program:
        parser.add_argument("--outdir", help="assemble outdir", required=True)
        parser.add_argument("--sample", help="vdj sample name", required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument('--fastq_dir', required=True)
    parser.add_argument('--mode', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)
    parser.add_argument('--thread', help='thread', default=20)

