import subprocess


import celescope.tools.utils as utils
from celescope.tools.step import s_common
from celescope.tools.mkref import parse_genomeDir

class StarMixin():
    """
    Mixin class for STAR
    """
    def __init__(self, args):
        self.fq = args.fq
        self.genomeDir = args.genomeDir
        self.out_unmapped = args.out_unmapped
        self.debug = args.debug
        self.outFilterMatchNmin = int(args.outFilterMatchNmin)
        self.multi_max = int(args.outFilterMultimapNmax)
        self.STAR_param = args.STAR_param
        self.consensus_fq = args.consensus_fq

        # parse
        self.genome = parse_genomeDir(self.genomeDir)

        # out 
        self.outPrefix = f'{self.outdir}/{self.sample}_'
        self.STAR_map_log = f'{self.outdir}/{self.sample}_Log.final.out'
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
            '--outSAMtype', 'BAM', 'Unsorted', # controls sort by Coordinate or not
            '--outFilterMatchNmin', str(self.outFilterMatchNmin)
        ]
        if self.out_unmapped:
            cmd += ['--outReadsUnmapped', 'Fastx']
        if self.fq[-3:] == ".gz":
            cmd += ['--readFilesCommand', 'zcat']
        cmd = ' '.join(cmd)
        if self.STAR_param:
            cmd += (" " + self.STAR_param)
        StarMixin.STAR.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run_star(self):
        self.STAR()
        self.sort_bam()
        self.index_bam()

    @utils.add_log
    def sort_bam(self):
        cmd = (
            f'samtools sort {self.unsort_STAR_bam} '
            f'-o {self.STAR_bam} '
            f'--threads {self.thread} '
        )
        StarMixin.sort_bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def index_bam(self):
        cmd = f"samtools index {self.STAR_bam}"
        StarMixin.index_bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)


def get_opts_star_mixin(parser, sub_program):
    parser.add_argument('--outFilterMatchNmin', help='STAR outFilterMatchNmin', default=0)
    parser.add_argument('--out_unmapped', help='out_unmapped', action='store_true')
    parser.add_argument('--genomeDir', help='genome directory')
    parser.add_argument('--STAR_param', help='STAR parameters', default="")
    parser.add_argument('--outFilterMultimapNmax', help='STAR outFilterMultimapNmax', default=1)
    parser.add_argument('--starMem', help='starMem', default=30)
    if sub_program:
        parser.add_argument('--fq', required=True)
        parser.add_argument("--consensus_fq", action='store_true', help="input fastq is umi consensus")
        parser = s_common(parser)
