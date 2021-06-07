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


class Go_assemble(Step):
    """
    Features

    - Assemble TCR/BCR full length by tracer.
    - Summary mapping rate.

    Output

    - `04.go_assemble/tracer` or `04.go_assemble/bracer` Tracer output directory.
    - `04.go_assemble/stat.txt` Recording mapping rate.
    """
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        self.species = args.species
        self.type = args.type
        self.thread = int(args.thread)
        self.fastq_dir = args.fastq_dir


    def run_tracer(self):

        fqs = [join(self.fastq_dir, f) for f in listdir(self.fastq_dir) if isfile(join(self.fastq_dir, f))]
        outdirs = [self.outdir] * len(fqs)
        species = [self.species] * len(fqs)
        if not os.path.exists(f'{self.outdir}/tracer'):
            os.makedirs(f'{self.outdir}/tracer')

        all_res = []
        with ProcessPoolExecutor(self.thread) as pool:
            for res in pool.map(tracer, fqs, outdirs, species):
                all_res.append(res)

        tracer_summarise(self.outdir)


    def run_bracer(self):
        fqs = [join(self.fastq_dir, f) for f in listdir(self.fastq_dir) if isfile(join(self.fastq_dir, f))]
        outdirs = [self.outdir] * len(fqs)
        species = [self.species] * len(fqs)
        if not os.path.exists(f'{self.outdir}/bracer'):
            os.makedirs(f'{self.outdir}/bracer')

        all_res = []
        with ProcessPoolExecutor(self.thread) as pool:
            for res in pool.map(bracer, fqs, outdirs, species):
                all_res.append(res)

        bracer_summarise(self.outdir)


################def get_reads_count(fq):
#    with pysam.FastxFile(fq) as fh:
#        count = 0
#        for entry in fh:
#            count += 1
#    return count

    def go_assemble_summary(self):

        count_file = f'{self.outdir}/../03.split_fastq/{self.sample}_count.txt'
        UMIs = pd.read_csv(count_file, sep='\t')

        all_UMIs = UMIs['UMIs_count'].tolist()
        medians = int(np.median(all_UMIs))
        all_UMIs = sum(all_UMIs)
        
        stat_file = self.outdir + '/stat.txt'

        go_assemble_summary = []

        if type == 'TCR':
            TRAs = glob.glob(f'{self.outdir}/tracer/*/aligned_reads/*_TCR_A.fastq')
            TRBs = glob.glob(f'{self.outdir}/tracer/*/aligned_reads/*_TCR_B.fastq')
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

            with open(f'{self.outdir}/tmp.txt', 'w') as f:
                f.write(f'Madian UMIs per cell:{medians}\n')
                f.write(f'Median TRA UMIs per cell:{medianA}\n')
                f.write(f'Median TRB UMIs per cell:{medianB}\n')

        elif type == 'BCR':
            IGHs = glob.glob(f'{self.outdir}/bracer/*/aligned_reads/*_BCR_H.fastq')
            IGKs = glob.glob(f'{self.outdir}/bracer/*/aligned_reads/*_BCR_K.fastq')
            IGLs = glob.glob(f'{self.outdir}/bracer/*/aligned_reads/*_BCR_L.fastq')

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

            with open(f'{self.outdir}/tmp.txt', 'w') as f:
                f.write(f'Median UMIs per cell:{medians}\n')
                f.write(f'Median IGH UMIs per Cell:{medianH}\n')
                f.write(f'Median IGK UMIs per Cell:{medianK}\n') 
                f.write(f'Median IGL UMIs per Cell:{medianL}\n')
                
        df = pd.DataFrame(go_assemble_summary, columns=['item', 'count', 'total_count'])

        utils.gen_stat(df, stat_file)

        self.clean_up()


    @utils.add_log
    def run(self):
        if self.type == 'TCR':
            self.run_tracer()
        elif self.type == 'BCR':
            self.run_bracer()
        self.go_assemble_summary()


@utils.add_log
def go_assemble(args):
    step_name = 'go_assemble'
    go_assemble_obj = Go_assemble(args, step_name)
    go_assemble_obj.run()
    

def get_opts_go_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fastq_dir', required=True)
    parser.add_argument('--type', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)

