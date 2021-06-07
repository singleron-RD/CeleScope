import argparse
import os
from os import listdir
from os.path import isfile, join
from concurrent.futures import ProcessPoolExecutor
from celescope.tools import utils
from celescope.tools.utils import *
import datetime
import glob
import pysam
import numpy as np
from celescope.tools.Step import Step, s_common


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
        f'--no_trimming '
        f'-r '
        f'--species {species} '
        f'-c {BRACER_CONF} '
        f'{prefix} '
        f'{outdir}/bracer '
        f'{fq} '
    )
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
        f'-r '
        f'--species {species} '
        f'-c {CONF_PATH} '
        f'{fq} '
        f'{prefix} '
        f'{outdir}/tracer '
    )
    os.system(cmd)


@utils.add_log
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


@utils.add_log
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


################def get_reads_count(fq):
#    with pysam.FastxFile(fq) as fh:
#        count = 0
#        for entry in fh:
#            count += 1
#    return count


def get_umi_count(fq):
    umis = []
    with pysam.FastxFile(fq) as fh:
        for entry in fh:
            attr = entry.name.split('_')
            barcode = attr[0]
            umi = attr[1]
            umis.append(umi)
    res = len(set(umis))
    return res


def go_assemble_summary(outdir, type):

    total_fq = f'{outdir}/../03.split_fastq/reads_count.tsv'
    UMIs = pd.read_csv(total_fq, sep='\t')

    all_UMIs = UMIs['UMIs_count'].tolist()
    medians = int(np.median(all_UMIs))
    all_UMIs = sum(all_UMIs)
    
    stat_file = outdir + '/../04.go_assemble/stat.txt'

    go_assemble_summary = []

    if type == 'TCR':
        TRAs = glob.glob(f'{outdir}/tracer/*/aligned_reads/*_TCR_A.fastq')
        TRBs = glob.glob(f'{outdir}/tracer/*/aligned_reads/*_TCR_B.fastq')
        TRA_UMIs = [get_umi_count(fq) for fq in TRAs]
        TRB_UMIs = [get_umi_count(fq) for fq in TRBs]
        TRA_UMIs_count = sum(TRA_UMIs)
        medianA = int(np.median(TRA_UMIs))
        TRB_UMIs_count = sum(TRB_UMIs)
        medianB = int(np.median(TRB_UMIs))

        totals = TRA_UMIs_count + TRB_UMIs_count

        go_assemble_summary.append({
            'item': f'All UMIs mapped to TRA or TRB',
            'count': totals,
            'total_count': all_UMIs, 
        })

        go_assemble_summary.append({
            'item': f'UMIs mapped to TRA',
            'count': TRA_UMIs_count,
            'total_count': all_UMIs,
        })

        go_assemble_summary.append({
            'item': f'UMIs mapped to TRB',
            'count': TRB_UMIs_count,
            'total_count': all_UMIs,
        })

        with open(f'{outdir}/tmp.txt', 'w') as f:
            f.write(f'Madian UMIs per cell:{medians}\n')
            f.write(f'Median TRA UMIs per cell:{medianA}\n')
            f.write(f'Median TRB UMIs per cell:{medianB}\n')

    elif type == 'BCR':
        IGHs = glob.glob(f'{outdir}/bracer/*/aligned_reads/*_BCR_H.fastq')
        IGKs = glob.glob(f'{outdir}/bracer/*/aligned_reads/*_BCR_K.fastq')
        IGLs = glob.glob(f'{outdir}/bracer/*/aligned_reads/*_BCR_L.fastq')

        IGH_UMIs = [get_umi_count(fq) for fq in IGHs]
        IGK_UMIs = [get_umi_count(fq) for fq in IGKs]
        IGL_UMIs = [get_umi_count(fq) for fq in IGLs]

        IGH = sum(IGH_UMIs)
        medianH = np.median(IGH_UMIs)
        IGK = sum(IGK_UMIs)
        medianK = np.median(IGK_UMIs)
        IGL = sum(IGL_UMIs)
        medianL = np.median(IGL_UMIs)

        totals = IGH + IGK + IGL

        go_assemble_summary.append({
            'item': f'All UMIs mapped to IGH, IGL or IGK',
            'count': totals,
            'total_count': all_UMIs,            
        })

        go_assemble_summary.append({
            'item': f'UMIs mapped to IGH',
            'count': IGH,
            'total_count': all_UMIs,
        })

        go_assemble_summary.append({
            'item': f'UMIs mapped to IGK',
            'count': IGK,
            'total_count': all_UMIs,
        })

        go_assemble_summary.append({
            'item': f'UMIs mapped to IGL',
            'count': IGL,
            'total_count': all_UMIs,
        })

        with open(f'{outdir}/tmp.txt', 'w') as f:
            f.write(f'Median UMIs per cell:{medians}\n')
            f.write(f'Median IGH UMIs per Cell:{medianH}\n')
            f.write(f'Median IGK UMIs per Cell:{medianK}\n') 
            f.write(f'Median IGL UMIs per Cell:{medianL}\n')
            
    df = pd.DataFrame(go_assemble_summary, columns=['item', 'count', 'total_count'])

    utils.gen_stat(df, stat_file)


def go_assemble(args):
    step_name = 'go_assemble'
    step = Step(args, step_name)
    thread = int(args.thread)
    fastq_dir = args.fastq_dir
    outdir = args.outdir
    species = args.species
    
    type = args.type
    if type == 'TCR':
        run_tracer(outdir, fastq_dir, species, thread)
    elif type == 'BCR':
        run_bracer(outdir, fastq_dir, species, thread)

    go_assemble_summary(outdir, type)

    step.clean_up()

def get_opts_go_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fastq_dir', required=True)
    parser.add_argument('--type', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)

