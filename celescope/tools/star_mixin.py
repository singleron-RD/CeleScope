import re
import subprocess
import os

from celescope.tools import utils
from celescope.tools.mkref import Mkref
from celescope.tools.step import s_common
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step


@utils.add_log
def get_star_cmd(args, input_file, output_prefix):
    """
    output sam format to improve speed
    """
    cmd = (
        f'STAR '
        f'--runThreadN {args.thread} '
        f'--genomeDir {args.genomeDir} '
        f'--outFilterMultimapNmax {args.outFilterMultimapNmax} '
        f'--outSAMtype BAM Unsorted '
        f'--outFilterMatchNmin {args.outFilterMatchNmin} '
        f'{args.STAR_param} '
        f'--readFilesIn {input_file} '
        f'--outFileNamePrefix {output_prefix}_ '
    )
    get_star_cmd.logger.info(cmd)
    return cmd


class Star_mixin(Step):
    """
    Mixin class for STAR
    """

    def __init__(self, args, add_prefix=None, display_title=None):
        super().__init__(args, display_title)

        # parse
        self.genome = Mkref.parse_genomeDir(args.genomeDir)
        self.stat_prefix = 'Reads'
        if getattr(args, "consensus_fq", False):
            self.stat_prefix = 'UMIs'

        # out
        if add_prefix:
            self.out_prefix += f'_{add_prefix}'
        self.STAR_map_log = f'{self.out_prefix}_Log.final.out'
        self.unsort_STAR_bam = f'{self.out_prefix}_Aligned.out.bam'
        self.STAR_bam = f'{self.out_prefix}_Aligned.sortedByCoord.out.bam'

    @utils.add_log
    def STAR(self):
        input_file = self.args.fq
        output_prefix = self.out_prefix
        cmd = get_star_cmd(self.args, input_file, output_prefix)

        if self.args.out_unmapped:
            cmd += '--outReadsUnmapped Fastx'
        if self.args.fq[-3:] == ".gz":
            cmd += '--readFilesCommand zcat'
        if self.args.STAR_param:
            cmd += (" " + self.args.STAR_param)
        self.STAR.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.STAR()
        self.add_star_metrics()
        self.sort_bam()
        self.index_bam()
        self.remove_unsort_bam()

    @utils.add_log
    def sort_bam(self):
        utils.sort_bam(
            self.unsort_STAR_bam,
            self.STAR_bam,
            threads=self.thread,
        )

    @utils.add_log
    def index_bam(self):
        utils.index_bam(self.STAR_bam)

    @utils.add_log
    def remove_unsort_bam(self):
        os.remove(self.unsort_STAR_bam)

    def add_star_metrics(self):
        """
        step metrics
        """

        with open(self.STAR_map_log, 'r') as map_log:
            # number amd percent
            unique_reads_list = []
            multi_reads_list = []
            total_reads = 0
            for line in map_log:
                if line.strip() == '':
                    continue
                if re.search(r'Uniquely mapped reads', line):
                    unique_reads_list.append(line.strip().split()[-1])
                if re.search(r'of reads mapped to too many loci', line):
                    multi_reads_list.append(line.strip().split()[-1])
                if re.search(r'Number of input reads', line):
                    total_reads = int(line.strip().split()[-1])

        unique_reads = int(unique_reads_list[0])
        multi_reads = int(multi_reads_list[0])

        self.add_metric(
            name='Genome',
            value=self.genome['genome_name'],
        )
        self.add_metric(
            name=f'Uniquely Mapped {self.stat_prefix}',
            value=unique_reads,
            total=total_reads,
            help_info='reads that mapped uniquely to the genome'
        )
        self.add_metric(
            name=f'Multi-Mapped {self.stat_prefix}',
            value=multi_reads,
            total=total_reads,
            help_info='reads that mapped to multiple locations in the genome'
        )


def get_opts_star_mixin(parser, sub_program):
    parser.add_argument(
        '--genomeDir',
        help=HELP_DICT['genomeDir'],
    )
    parser.add_argument(
        '--outFilterMatchNmin',
        help="""Alignment will be output only if the number of matched bases 
is higher than or equal to this value.""",
        default=50,
    )
    parser.add_argument(
        '--out_unmapped',
        help='Output unmapped reads.',
        action='store_true'
    )
    parser.add_argument('--STAR_param', help=HELP_DICT['additional_param'], default="")
    parser.add_argument(
        '--outFilterMultimapNmax',
        help='Default `1`. How many places are allowed to match a read at most.',
        default=1
    )
    parser.add_argument(
        '--starMem',
        help='Default `30`. Maximum memory that STAR can use.',
        default=30
    )
    if sub_program:
        parser.add_argument('--fq', help="Required. R2 fastq file.", required=True)
        parser.add_argument("--consensus_fq", action='store_true',
                            help="A indicator that the input fastq has been consensused.")
        parser = s_common(parser)
