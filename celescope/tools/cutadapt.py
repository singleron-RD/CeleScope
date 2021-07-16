import re
import subprocess
from itertools import islice

import pandas as pd
import pysam

from celescope.tools.step import Step, s_common
import celescope.tools.utils as utils

ADAPTER = ['polyT=A{18}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC']


class Cutadapt(Step):
    """
    Features
    - Trim adapters in R2 reads with cutadapt. Default adapters includes:
        - polyT=A{18}, 18 A bases. 
        - p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Illumina p5 adapter.

    Output
    - `cutadapt.log` Cutadapt output log file.
    - `{sample}_clean_2.fq.gz` R2 reads file without adapters.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # set
        self.adapter_args = self.read_adapter_fasta(args.adapter_fasta)
        self.adapter_args += ADAPTER

        # out files
        if args.gzip:
            suffix = ".gz"
        else:
            suffix = ""
        self.out_fq2 = f'{self.outdir}/{self.sample}_clean_2.fq{suffix}'
        self.cutadapt_log_file = f'{self.outdir}/cutadapt.log'

    @staticmethod
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

    def format_and_write_stat(self, cutadapt_log):
        # Total reads processed:...Total written (filtered):
        start_line = 9
        end_line = 16
        useful_content = islice(cutadapt_log.split('\n'), start_line, end_line)
        p_list = []
        for line_index, line in enumerate(useful_content):
            if line.strip() == '':
                continue
            line = re.sub(r'\s{2,}', r'', line)
            line = re.sub(r' bp', r'', line)
            line = re.sub(r'(?<=\d)\s+\(', r'(', line)
            line = line.strip()
            attr = line.split(":")
            if len(attr) < 2:
                raise Exception(
                    'May caused by cutadapt version != 1.1.7. '
                    f'Check version in {self.outdir}/cutadapt.log\n'
                    f'cutadapt log error at line {line_index + start_line}:\n'
                    f'{line}'
                )
            p_list.append({"item": attr[0], "value": attr[1]})
        p_df = pd.DataFrame(p_list)
        p_df.iloc[0, 0] = 'Reads with Adapters'
        p_df.iloc[1, 0] = 'Reads too Short'
        p_df.iloc[2, 0] = 'Reads Written'
        p_df.iloc[3, 0] = 'Base Pairs Processed'
        p_df.iloc[4, 0] = 'Base Pairs Quality-Trimmed'
        p_df.iloc[5, 0] = 'Base Pairs Written'
        p_df.to_csv(self.stat_file, sep=':', index=False, header=None)

    @utils.add_log
    def run(self):
        adapter_args_str = " ".join(['-a ' + adapter for adapter in self.adapter_args])
        cmd = (
            'cutadapt '
            f'{adapter_args_str} '
            f'-n {len(self.adapter_args)} '
            f'-j {self.thread} '
            f'-m {self.args.minimum_length} '
            f'--nextseq-trim={self.args.nextseq_trim} '
            f'--overlap {self.args.overlap} '
            f'-l {self.args.insert} '
            f'-o {self.out_fq2} '
            f'{self.args.fq} '
        )
        Cutadapt.run.logger.info(cmd)
        # need encoding argument to return str
        results = subprocess.run(
            cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE,
            encoding='utf-8', check=True, shell=True
        )
        cutadapt_log = results.stdout
        with open(self.cutadapt_log_file, 'w') as f:
            f.write(cutadapt_log)
        self.format_and_write_stat(cutadapt_log)
        self.clean_up()


@utils.add_log
def cutadapt(args):

    step_name = "cutadapt"
    cutadapt_obj = Cutadapt(args, step_name)
    cutadapt_obj.run()


def get_opts_cutadapt(parser, sub_program):
    parser.add_argument('--adapter_fasta', help='Addtional adapter fasta file.')
    parser.add_argument(
        '--minimum_length',
        help='Default `20`. Discard processed reads that are shorter than LENGTH.',
        default=20
    )
    parser.add_argument(
        '--nextseq_trim',
        help="""Default `20`. Quality trimming of reads using two-color chemistry (NextSeq). 
Some Illumina instruments use a two-color chemistry to encode the four bases. 
This includes the NextSeq and the NovaSeq. 
In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
However, dark cycles also occur when sequencing “falls off” the end of the fragment.
The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.""",
        default=20,
    )
    parser.add_argument(
        '--overlap',
        help="""Default `10`. Since Cutadapt allows partial matches between the read and the adapter sequence,
short matches can occur by chance, leading to erroneously trimmed bases. 
For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter. 
To reduce the number of falsely trimmed bases, the alignment algorithm requires that 
at least {overlap} bases match between adapter and read. """,
        default=10
    )
    parser.add_argument(
        '--insert',
        help="Default `150`. Read2 insert length.",
        default=150
    )
    if sub_program:
        parser.add_argument('--fq', help='Required. R2 reads from step Barcode.', required=True)
        parser.add_argument('--gzip', help="Output gzipped fastq", action='store_true')
        parser = s_common(parser)
    return parser
