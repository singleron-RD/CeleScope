import re
import subprocess

import celescope.tools.utils as utils
from celescope.tools.mkref import Mkref
from celescope.tools.step import s_common
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step


class Star_mixin(Step):
    """
    Mixin class for STAR
    """

    def __init__(self, args, add_prefix=None, display_title=None):
        super().__init__(args, display_title)

        self.fq = args.fq
        self.genomeDir = args.genomeDir
        self.out_unmapped = args.out_unmapped
        self.debug = args.debug
        self.outFilterMatchNmin = int(args.outFilterMatchNmin)
        self.multi_max = int(args.outFilterMultimapNmax)
        self.STAR_param = args.STAR_param
        self.consensus_fq = args.consensus_fq

        # parse
        self.genome = Mkref.parse_genomeDir(self.genomeDir)
        self.stat_prefix = 'Reads'
        if self.consensus_fq:
            self.stat_prefix = 'UMIs'

        # out
        self.outPrefix = f'{self.outdir}/{self.sample}_'
        if add_prefix:
            self.outPrefix += add_prefix + '_'
        self.STAR_map_log = f'{self.outPrefix}Log.final.out'
        self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
        self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'

    @utils.add_log
    def STAR(self):
        cmd = [
            'STAR',
            '--runThreadN', str(self.thread),
            '--genomeDir', self.genomeDir,
            '--readFilesIn', self.fq,
            '--outFilterMultimapNmax', str(self.multi_max),
            '--outFileNamePrefix', self.outPrefix,
            '--outSAMtype', 'BAM', 'Unsorted',  # controls sort by Coordinate or not
            '--outFilterMatchNmin', str(self.outFilterMatchNmin)
        ]
        if self.out_unmapped:
            cmd += ['--outReadsUnmapped', 'Fastx']
        if self.fq[-3:] == ".gz":
            cmd += ['--readFilesCommand', 'zcat']
        cmd = ' '.join(cmd)
        if self.STAR_param:
            cmd += (" " + self.STAR_param)
        self.STAR.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.STAR()
        self.get_star_metrics()
        self.sort_bam()
        self.index_bam()

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

    def get_star_metrics(self):
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
        help="""Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.""",
        default=0
    )
    parser.add_argument(
        '--out_unmapped',
        help='Output unmapped reads.',
        action='store_true'
    )
    parser.add_argument('--STAR_param', help='Other STAR parameters.', default="")
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
