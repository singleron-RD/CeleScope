import re
import subprocess

import pysam

from celescope.tools.step import Step, s_common
from celescope.tools import utils

ADAPTER = ['polyT=A{18}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC']

LOG_METRICS_TITLE = (
    'Total reads processed',
    'Reads with adapters',
    'Reads that were too short',
    'Reads written (passing filters)',
    'Total basepairs processed',
    'Quality-trimmed',
    'Total written (filtered)',
)

def read_cutadapt_log(cutadapt_log):
    """
    Returns:
        dict: {
            'Total reads processed': int, 'Reads with adapters': int, 'Reads that were too short': int, 
            'Reads written (passing filters)': int, 'Total basepairs processed': int, 
            'Quality-trimmed': int, 'Total written (filtered)': int
        }
    """
    metrics_dict = {}
    remove_strs = [r',', r' bp', r'\(.*\)']

    for _line_index, line in enumerate(cutadapt_log.split('\n')):
        line = line.strip()
        if not line:
            continue
        attr = line.split(":")
        if attr[0] in LOG_METRICS_TITLE:
            title, number = attr[0], attr[1]                
            number = attr[1].strip()
            for remove_str in remove_strs:
                number = re.sub(pattern=remove_str, repl='', string=number)
            metrics_dict[title] = int(number)

    return metrics_dict

class Cutadapt(Step):
    """
    ## Features
    - Trim adapters in R2 reads with cutadapt. Default adapters includes:
        - polyT=A{18}, 18 A bases. 
        - p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Illumina p5 adapter.
    ## Output
    - `cutadapt.log` Cutadapt output log file.
    - `{sample}_clean_2.fq.gz` R2 reads file without adapters.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

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



    def add_cutadapt_metrics(self, cutadapt_log):
        # Total reads processed:...Total written (filtered):
        metrics_dict = read_cutadapt_log(cutadapt_log)

        total_reads = metrics_dict['Total reads processed']
        reads_with_adapters = metrics_dict['Reads with adapters']
        reads_too_short = metrics_dict['Reads that were too short']
        reads_written = metrics_dict['Reads written (passing filters)']
        total_base_pairs = metrics_dict['Total basepairs processed']
        quality_trimmed = metrics_dict['Quality-trimmed']
        base_pairs_written = metrics_dict['Total written (filtered)']

        self.add_metric(
            name='Reads with Adapters',
            value=reads_with_adapters,
            total=total_reads,
            help_info='R2 reads with poly A(read-through adpaters) or sequencing adapters'
        )
        self.add_metric(
            name='Reads too Short',
            value=reads_too_short,
            total=total_reads,
            help_info='reads with read length less than 20bp after trimming'
        )
        self.add_metric(
            name='Reads Written',
            value=reads_written,
            total=total_reads,
            help_info='reads pass filtering'
        )
        self.add_metric(
            name='Base Pairs Processed',
            value=total_base_pairs,
            help_info='total processed base pairs'
        )
        self.add_metric(
            name='Base Pairs Quality-Trimmed',
            value=quality_trimmed,
            total=total_base_pairs,
            help_info='bases pairs removed from the end of the read whose quality is smaller than the given threshold'
        )
        self.add_metric(
            name='Base Pairs Written',
            value=base_pairs_written,
            total=total_base_pairs,
            help_info='base pairs pass filtering'
        )

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
            f'{self.args.cutadapt_param} '
            f'{self.args.fq} '
        )
        self.run.logger.info(cmd)
        # need encoding argument to return str
        results = subprocess.run(
            cmd, stdout=subprocess.PIPE,
            encoding='utf-8', check=True, shell=True
        )
        cutadapt_log = results.stdout
        with open(self.cutadapt_log_file, 'w') as f:
            f.write(cutadapt_log)
        self.add_cutadapt_metrics(cutadapt_log)


@utils.add_log
def cutadapt(args):
    with Cutadapt(args, display_title="Trimming") as runner:
        runner.run()


def get_opts_cutadapt(parser, sub_program):
    parser.add_argument('--gzip', help="Output gzipped fastq files.", action='store_true')
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
    parser.add_argument('--cutadapt_param', help='Other cutadapt parameters. For example, --cutadapt_param "-g AAA" ', default="")
    if sub_program:
        parser.add_argument('--fq', help='Required. R2 reads from step Barcode.', required=True)
        parser = s_common(parser)
    return parser
