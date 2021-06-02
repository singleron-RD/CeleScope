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


def get_assemble_stat(outdir, type):

    total_fq = f'{outdir}/../03.split_fastq/reads_count.tsv'
    UMIs = pd.DataFrame(total_fq, sep='\t')

    all_UMIs = UMIs['UMIs_count'].sum()
    stat_file = outdir + '/../04.go_assemble/stat.txt'

    if type == 'TCR':
        TRAs = glob.glob(f'{outdir}/tracer/*/aligned_reads/*_TCR_A.fastq')
        TRBs = glob.glob(f'{outdir}/tracer/*/aligned_reads/*_TCR_B.fastq')
        TRA_UMIs = [get_umi_count(fq) for fq in TRAs]
        TRB_UMIs = [get_umi_count(fq) for fq in TRBs]
        TRA_UMIs_count = sum(TRA_UMIs)
        TRA_ = format(TRA_UMIs_count, ',')
        TRB_UMIs_count = sum(TRB_UMIs)
        TRB_ = format(TRB_UMIs_count, ',')

        TRA_mapping = TRA_UMIs_count/all_UMIs
        TRA_mapping = round(TRA_mapping, 4)
        TRA_mapping = f'{TRA_}({TRA_mapping})'

        TRB_mapping = TRB_UMIs_count/all_UMIs
        TRB_mapping = round(TRB_mapping, 4)
        TRB_mapping = f'{TRB_}({TRB_mapping})'

        total_counts = TRA_UMIs_count + TRB_UMIs_count
        total_ = format(total_counts, ',')
        total_mapping = (total_counts)/all_UMIs
        total_mapping = round(total_mapping, 4)
        total_mapping = f'{total_}({total_mapping})'

        stat_text = pd.DataFrame({
            'item': ['UMIs mapped to TRA or TRB', 'UMIs mapped to TRA', 'UMIs mapped to TRB'], 'count': [total_mapping, TRA_mapping, TRB_mapping]
        }, columns=['item', 'count'])

        stat_text.to_csv(stat_file, sep=':', header=None, index=False)

    elif type == 'BCR':
        IGHs = glob.glob(f'{outdir}/bracer/*/aligned_reads/*_BCR_H.fastq')
        IGKs = glob.glob(f'{outdir}/bracer/*/aligned_reads/*_BCR_K.fastq')
        IGLs = glob.glob(f'{outdir}/bracer/*/aligned_reads/*_BCR_L.fastq')

        IGH_UMIs = [get_umi_count(fq) for fq in IGHs]
        IGK_UMIs = [get_umi_count(fq) for fq in IGKs]
        IGL_UMIs = [get_umi_count(fq) for fq in IGLs]


        IGH = sum(IGH_UMIs)
        IGH_ = format(IGH, ',')
        IGK = sum(IGK_UMIs)
        IGK_ = format(IGK, ',')
        IGL = sum(IGL_UMIs)
        IGL_ = format(IGL, ',')

        IGH_mapping = IGH/all_UMIs
        IGH_mapping = round(IGH_mapping, 4)
        IGH_mapping = f'{IGH_}({IGH_mapping})'

        IGK_mapping = IGK/all_UMIs
        IGK_mapping = round(IGK_mapping, 4)
        IGK_mapping = f'{IGK_}({IGK_mapping})'

        IGL_mapping = IGL/all_UMIs
        IGL_mapping = round(IGL_mapping, 4)
        IGL_mapping = f'{IGL_}({IGL_mapping})'

        total_counts = IGH + IGK + IGL
        total_ = format(total_counts, ',')

        total_mapping = (total_counts)/all_UMIs
        total_mapping = round(total_mapping, 4)
        total_mapping = f'{total_}({total_mapping})'

        stat_text = pd.DataFrame({
            'item': ['UMIs mapped to IGH, IGK or IGL', 'UMIs mapped to IGH', 'UMIs mapped to IGK', 'UMIs mapped to IGL'], 'count': [total_mapping, IGH_mapping, IGK_mapping, IGL_mapping]
        })

        stat_text.to_csv(stat_file, sep=':', header=None, index=False)






def go_assemble(args):
    thread = int(args.thread)
    fastq_dir = args.fastq_dir
    outdir = args.outdir
    species = args.species
    
    type = args.type
    if type == 'TCR':
        run_tracer(outdir, fastq_dir, species, thread)
    elif type == 'BCR':
        run_bracer(outdir, fastq_dir, species, thread)


def get_opts_go_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fastq_dir', required=True)
    parser.add_argument('--type', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)

