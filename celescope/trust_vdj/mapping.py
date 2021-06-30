import re
import subprocess
from collections import defaultdict
import celescope.tools.utils as utils
from celescope.tools.mkref import parse_genomeDir
from celescope.tools.step import s_common, Step
import pysam
import pandas as pd


MAPPING_INDEX = '/SGRNJ03/randd/zhouxin/software/TRUST4/index'


class Mapping(Step):
    """
    Mixin class for STAR
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        self.fq = args.fq
        self.out_unmapped = args.out_unmapped
        self.debug = args.debug
        self.outFilterMatchNmin = int(args.outFilterMatchNmin)
        self.multi_max = int(args.outFilterMultimapNmax)
        self.STAR_param = args.STAR_param
        self.consensus_fq = args.consensus_fq
        self.Seqtype = args.Seqtype
        self.species = args.species

        # parse

        self.stat_prefix = 'Reads'
        if self.consensus_fq:
            self.stat_prefix = 'UMIs'

        # out
        self.outPrefix = f'{self.outdir}/{self.sample}_'
        self.STAR_map_log = f'{self.outPrefix}Log.final.out'
        self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
        self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'

    @utils.add_log
    def STAR(self):
        genome = f'{MAPPING_INDEX}/{self.species}'
        cmd = [
            'STAR',
            '--runThreadN', str(self.thread),
            '--genomeDir', genome,
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
        Mapping.STAR.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
    

    def run_star(self):
        self.STAR()
        self.sort_bam()
        self.index_bam()
        self.get_star_metrics()


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
    def get_star_metrics(self):
        """
        step metrics
        """

        mapping_summary = []
        with open(self.STAR_map_log, 'r') as map_log:
            # number amd percent
            unique_reads_list = []
            multi_reads_list = []
            for line in map_log:
                if line.strip() == '':
                    continue
                if re.search(r'Uniquely mapped reads', line):
                    unique_reads_list.append(line.strip().split()[-1])
                if re.search(r'of reads mapped to too many loci', line):
                    multi_reads_list.append(line.strip().split()[-1])
                if 'Number of input reads' in line:
                    all_reads_list = re.findall(r"\d+", line)
                    all_reads = int(all_reads_list[0])

        with pysam.AlignmentFile(self.STAR_bam) as bam:
            dic = defaultdict(list)
            if self.Seqtype == 'TCR':
                loci = ['TRA', 'TRB']
                prefix = 'TR'
            elif self.Seqtype == 'BCR':
                prefix = 'IG'
                loci = ['IGH', 'IGL', 'IGK']
            for read in bam:
                if prefix in read.reference_name:
                    dic[prefix].append(read)
                for l in loci:
                    if l in read.reference_name:
                        dic[l].append(read)

            total_count = len(dic[prefix])

            mapping_summary.append({
                'item': f'Reads mapped to any {self.Seqtype} V(D)J genes',
                'count': total_count,
                'total_count': all_reads,
                }
            )
            for l in loci:
                tmp_count = len(dic[l])
                tmp_name = f'Reads mapped to {l}'
                mapping_summary.append({
                    'item': tmp_name,
                    'count': tmp_count,
                    'total_count': all_reads
                })

            sum_df = pd.DataFrame(mapping_summary, columns=['item', 'count', 'total_count'])

            stat_file = self.outdir + '/stat.txt'
            utils.gen_stat(sum_df, stat_file)

            self.clean_up


def mapping(args):
    step_name = 'mapping'
    mapping_obj = Mapping(args, step_name)
    mapping_obj.run_star()

def get_opts_mapping(parser, sub_program):
    parser.add_argument(
        '--species',
        help='Required. Species name.', required=True, 
        choices=['Hsap', 'Mmus']
    )
    parser.add_argument(
        '--outFilterMatchNmin',
        help="""Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.""",
        default=0
    )
    parser.add_argument(
        '--out_unmapped',
        help='Output unmapped reads',
        action='store_true'
    )
    parser.add_argument('--STAR_param', help='Other STAR parameters', default="")
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
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        parser.add_argument('--fq', help="Required. R2 fastq file.", required=True)
        parser.add_argument("--consensus_fq", action='store_true', help="Input fastq has been consensused")
        parser = s_common(parser)


