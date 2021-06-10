import re
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
from concurrent.futures import ProcessPoolExecutor
from celescope.tools import utils
from celescope.tools.Step import Step, s_common


TRACER_PATH = '/SGRNJ03/randd/zhouxin/software/tracer/tracer'
CONF_PATH = '/SGRNJ03/randd/zhouxin/software/tracer/tracer.conf'
BRACER_PATH = '/SGRNJ03/randd/zhouxin/software/bracer/bracer'
BRACER_CONDA = 'bracer'
BRACER_CONF = '/SGRNJ03/randd/zhouxin/software/bracer/bracer.conf'


@utils.add_log
def assemble_summary(outdir, Seqtype, sample, species):
    
    stat_file = outdir + '/stat.txt'

    go_assemble_summary = []

    clean_fq = f'{outdir}/../02.cutadapt/{sample}_clean_2.fq'

    count_file = f'{outdir}/../03.split_fastq/{sample}_count.txt'

    count_df = pd.read_csv(count_file, sep='\t')

    total_count = count_df['readcount'].sum()

    if Seqtype == 'TCR':
        loci = ['A', 'B']

        total_mapped = 0

        for locus in loci:
            cmd = (
                f'source activate {BRACER_CONDA}; '
                f'bowtie2 -p 5 -k 1 --np 0 --rdg 1,1 --rfg 1,1 '
                f'-x /SGRNJ03/randd/zhouxin/software/tracer/resources/{species}/combinatorial_recombinomes/TCR_{locus} '
                f'-U {clean_fq} '
                f'-S {outdir}/TR{locus}.sam > {outdir}/log 2>&1'
            )
            os.system(cmd)
            with open(f'{outdir}/log') as fh:
                for line in fh:
                    if 'aligned exactly 1 time' in line:
                        res = re.findall(r"\d+", line)
                        item = f'Reads mapped to TR{locus}'
                        count = int(res[0])
                        total_mapped += count
                        go_assemble_summary.append({
                            'item': item,
                            'count': count,
                            'total_count': total_count,
                        })

            os.system(f'rm {outdir}/TR{locus}.sam')

        go_assemble_summary.insert(0, {
            'item': 'All reads Mapped to TRA and TRB',
            'count': total_mapped,
            'total_count': total_count
        })

        os.system(f'rm {outdir}/log')

    elif Seqtype == 'BCR':
        loci = ['H', 'L', 'K']

        total_mapped = 0

        for locus in loci:
            cmd = (
                f'source activate {BRACER_CONDA}; '
                f'bowtie2 -p 5 -k 1 --np 0 --rdg 1,1 --rfg 1,1 '
                f'-x /SGRNJ03/randd/zhouxin/software/bracer/resources/{species}/combinatorial_recombinomes/BCR_{locus} '
                f'-U {clean_fq} '
                f'-S {outdir}/BR{locus}.sam > {outdir}/log 2>&1'
            )
            os.system(cmd)
            with open(f'{outdir}/log') as fh:
                for line in fh:
                    if 'aligned exactly 1 time' in line:
                        res = re.findall(r"\d+", line)
                        item = f'Reads mapped to BR{locus}'
                        count = int(res[0])
                        total_mapped += count
                        go_assemble_summary.append({
                            'item': item,
                            'count': count,
                            'total_count': total_count,
                        })
            os.system(f'rm {outdir}/BR{locus}.sam')
        go_assemble_summary.insert(0, {
            'item': 'All reads Mapped to IGH, IGL and IGK',
            'count': total_mapped,
            'total_count': total_count
        })
        os.system(f'rm {outdir}/log')

    df = pd.DataFrame(go_assemble_summary, columns=['item', 'count', 'total_count'])

    utils.gen_stat(df, stat_file)


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
        self.Seqtype = args.Seqtype
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

        assemble_summary(self.outdir, self.Seqtype, self.sample, self.species)


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

        assemble_summary(self.outdir, self.Seqtype, self.sample, self.species)

    @utils.add_log
    def run(self):
        if self.Seqtype == 'TCR':
            self.run_tracer()
        elif self.Seqtype == 'BCR':
            self.run_bracer()

        self.clean_up()


@utils.add_log
def go_assemble(args):
    step_name = 'go_assemble'
    go_assemble_obj = Go_assemble(args, step_name)
    go_assemble_obj.run()
    

def get_opts_go_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fastq_dir', required=True)
    parser.add_argument('--Seqtype', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)

