import re
import subprocess
from collections import defaultdict
import celescope.tools.utils as utils
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
        self.name_sorted_bam = f'{self.outPrefix}name_sorted.bam'

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
        self.name_sort_bam()
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
    def name_sort_bam(self):
        cmd = (
            'samtools sort -n '
            f'-@ {self.thread} '
            f'-o {self.name_sorted_bam} '
            f'{self.STAR_bam} '
        )
        Mapping.name_sort_bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
            

    @utils.add_log
    def get_star_metrics(self):
        """
        step metrics
        """
        if self.Seqtype == 'TCR':
            chains = ['TRA', 'TRB']
        elif self.Seqtype == 'BCR':
            chains = ['IGH', 'IGL', 'IGK']

        mapping_summary = []
        with open(self.STAR_map_log, 'r') as map_log:
            # number amd percent
            for line in map_log:
                if line.strip() == '':
                    continue
                if 'Uniquely mapped reads number' in line:
                    unique_reads_list = re.findall(r"\d+", line)
                    unique_reads = int(unique_reads_list[0])
                if 'Number of reads mapped to multiple loci' in line:
                    multi_reads_list = re.findall(r"\d+", line)
                    multi_reads = int(multi_reads_list[0])
                if 'Number of reads mapped to too many loci' in line:
                    many_reads_list = re.findall(r"\d+", line)
                    many_reads = int(many_reads_list[0])
                if 'Number of input reads' in line:
                    all_reads_list = re.findall(r"\d+", line)
                    all_reads = int(all_reads_list[0])

        total_count = unique_reads+multi_reads+many_reads


        with pysam.AlignmentFile(self.name_sorted_bam) as bamfile:
            read_ref_dict = defaultdict(lambda: defaultdict(int))
            count_dic = defaultdict(int)
            for read in bamfile:
                name = read.query_name
                ref = read.reference_name
                for c in chains:
                    if ref.startswith(c):
                        read_ref_dict[name][c]+=1

            totals_ = 0
            for c in chains:
                totals_+=read_ref_dict[name][c]

            for c in chains:
                for name in read_ref_dict:
                    count = read_ref_dict[name][c]
                    if (count / totals_)*100 >= 80:
                        count_dic[c]+=1
                    else:
                        continue


            mapping_summary.append({
                'item': f'Reads mapped to any {self.Seqtype} V(D)J genes',
                'count': total_count,
                'total_count': all_reads,
                }
            )
            for c in chains:
                tmp_count = count_dic[c]
                tmp_name = f'Reads mapped to {c}'
                mapping_summary.append({
                    'item': tmp_name,
                    'count': tmp_count,
                    'total_count': all_reads
                })

            sum_df = pd.DataFrame(mapping_summary, columns=['item', 'count', 'total_count'])

            stat_file = self.outdir + '/stat.txt'
            utils.gen_stat(sum_df, stat_file)

            self.clean_up()


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


