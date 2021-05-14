import os
import re
import subprocess
import pandas as pd
import pysam

from itertools import islice
from celescope.tools.utils import add_log, s_common
from celescope.tools.step import Step


ADAPT = ['polyT=A{18}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC']


def format_stat(cutadapt_log):
    fh = open(cutadapt_log, 'r')
    stat_file = os.path.dirname(cutadapt_log) + '/stat.txt'
    # Total reads processed:...Total written (filtered):
    content = islice(fh, 9, 16)
    p_list = []
    for line in content:
        if line.strip() == '':
            continue
        line = re.sub(r'\s{2,}', r'', line)
        line = re.sub(r' bp', r'', line)
        line = re.sub(r'(?<=\d)\s+\(', r'(', line)
        line = line.strip()
        attr = line.split(":")
        p_list.append({"item": attr[0], "value": attr[1]})
    p_df = pd.DataFrame(p_list)
    p_df.iloc[0, 0] = 'Reads with Adapters'
    p_df.iloc[1, 0] = 'Reads too Short'
    p_df.iloc[2, 0] = 'Reads Written'
    p_df.iloc[3, 0] = 'Base Pairs Processed'
    p_df.iloc[4, 0] = 'Base Pairs Quality-Trimmed'
    p_df.iloc[5, 0] = 'Base Pairs Written'
    p_df.to_csv(stat_file, sep=':', index=False, header=None)

    fh.close()


@add_log
def read_adapter_fasta(adapter_fasta):
    '''
    return ['adapter1=AAA','adapter2=BBB']
    '''
    adapter_args = []
    if adapter_fasta and adapter_fasta != 'None':
        with pysam.FastxFile(adapter_fasta) as fh:
            for read in fh:
                adapter_args.append(f'{read.name}={read.sequence}')
    return adapter_args


@add_log
def cutadapt(args):

    step_name = "cutadapt"
    step = Step(args, step_name)

    adapter_args = read_adapter_fasta(args.adapter_fasta)
    adapter_args += ADAPT

    # run cutadapt
    adapt = []
    for a in adapter_args:
        adapt.append('-a')
        adapt.append(a)

    if args.gzip:
        suffix = ".gz"
    else:
        suffix = ""
    out_fq2 = f'{args.outdir}/{args.sample}_clean_2.fq{suffix}'
    cmd = ['cutadapt'] + adapt + ['-n',
                                  str(len(adapter_args)),
                                  '-j',
                                  str(args.thread),
                                  '-m',
                                  str(args.minimum_length),
                                  '--nextseq-trim=' + str(args.nextseq_trim),
                                  '--overlap',
                                  str(args.overlap),
                                  '-l',
                                  str(args.insert),
                                  '-o',
                                  out_fq2,
                                  args.fq]
    cutadapt.logger.info('%s' % (' '.join(cmd)))
    res = subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, check=True)
    with open(args.outdir + '/cutadapt.log', 'wb') as fh:
        fh.write(res.stdout)

    format_stat(args.outdir + '/cutadapt.log')

    step.clean_up()


def get_opts_cutadapt(parser, sub_program):
    parser.add_argument('--adapter_fasta', help='addtional adapter fasta file')
    parser.add_argument('--minimum_length',dest='minimum_length',help='minimum_length', default=20)
    parser.add_argument('--nextseq_trim', help='nextseq_trim', default=20)
    parser.add_argument('--overlap',help='minimum overlap length', default=10)
    parser.add_argument('--insert', help="read2 insert length", default=150)
    if sub_program:
        parser.add_argument('--fq', help='fq file', required=True)
        parser.add_argument('--gzip', help="output gzipped fastq", action='store_true')
        parser = s_common(parser)
    return parser


