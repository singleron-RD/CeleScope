import re
import subprocess
import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.utils import format_number, format_stat, glob_genomeDir
from celescope.tools.step import Step, s_common
from celescope.__init__ import ROOT_PATH


class Star(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        self.fq = args.fq
        self.genomeDir = args.genomeDir
        self.out_unmapped = args.out_unmapped
        self.debug = args.debug
        self.outFilterMatchNmin = int(args.outFilterMatchNmin)
        self.multi_max = int(args.outFilterMultimapNmax)
        self.STAR_param = args.STAR_param
        self.consensus_fq = args.consensus_fq

        if self.genomeDir and self.genomeDir != "None":
            self.STAR_index = self.genomeDir
        else:
            self.STAR_index = args.STAR_index
            self.refFlat = args.refFlat

        # out 
        self.outPrefix = f'{self.outdir}/{self.sample}_'
        self.STAR_map_log = f'{self.outdir}/{self.sample}_Log.final.out'
        self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
        self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'
        self.ribo_log = f'{self.outdir}/{self.sample}_ribo_log.txt'
        self.ribo_run_log = f'{self.outdir}/{self.sample}_ribo_run.log'
        self.picard_region_log = f'{self.outdir}/{self.sample}_region.log'
        self.plot = None
        self.stats = pd.Series()

    def format_stat(self):

        stat_prefix = 'Reads'
        if self.consensus_fq:
            stat_prefix = 'UMIs'

        fh1 = open(self.STAR_map_log, 'r')
        UNIQUE_READS = []
        MULTI_MAPPING_READS = []
        for line in fh1:
            if line.strip() == '':
                continue
            if re.search(r'Uniquely mapped reads', line):
                UNIQUE_READS.append(line.strip().split()[-1])
            if re.search(r'of reads mapped to too many loci', line):
                MULTI_MAPPING_READS.append(line.strip().split()[-1])
        fh1.close()

        fh2 = open(self.picard_region_log, 'r')
        region_dict = {}
        while True:
            line = fh2.readline()
            if not line:
                break
            if line.startswith('## METRICS CLASS'):
                header = fh2.readline().strip().split('\t')
                data = fh2.readline().strip().split('\t')
                region_dict = dict(zip(header, data))
                break
        fh2.close()
        
        Total = float(region_dict['PF_ALIGNED_BASES'])
        Exonic_Regions = int(region_dict['UTR_BASES']) + \
            int(region_dict['CODING_BASES'])
        Intronic_Regions = int(region_dict['INTRONIC_BASES'])
        Intergenic_Regions = int(region_dict['INTERGENIC_BASES'])

        region_dict['Exonic_Regions'] = "{}({:.2%})".format(
            format_number(Exonic_Regions), Exonic_Regions / Total)
        region_dict['Intronic_Regions'] = "{}({:.2%})".format(
            format_number(Intronic_Regions), Intronic_Regions / Total)
        region_dict['Intergenic_Regions'] = "{}({:.2%})".format(
            format_number(Intergenic_Regions), Intergenic_Regions / Total)

        self.stats = self.stats.append(pd.Series(
            f'{format_number(int(UNIQUE_READS[0]))}({UNIQUE_READS[1]})',
            index=[f'Uniquely Mapped {stat_prefix}']
        ))
        self.stats = self.stats.append(pd.Series(
            f'{format_number(int(MULTI_MAPPING_READS[0]))}({MULTI_MAPPING_READS[1]})',
            index=[f'Multi-Mapped {stat_prefix}']
        ))
        # ribo
        if self.debug:
            f = open(self.ribo_log, 'r')
            for line in f:
                if line.find('#Matched') != -1:
                    items = line.split()
                    Reads_Mapped_to_rRNA = int(items[1])
                if line.find('#Total') != -1:
                    items = line.split()
                    Reads_Total = int(items[1])

            self.stats = self.stats.append(pd.Series(
                format_stat(Reads_Mapped_to_rRNA, Reads_Total),
                index=[f'{stat_prefix} Mapped to rRNA']
            ))
            f.close()

        self.stats = self.stats.append(pd.Series(
            region_dict['Exonic_Regions'],
            index=['Base Pairs Mapped to Exonic Regions']
        ))
        self.stats = self.stats.append(pd.Series(
            region_dict['Intronic_Regions'],
            index=['Base Pairs Mapped to Intronic Regions']
        ))
        self.stats = self.stats.append(pd.Series(
            region_dict['Intergenic_Regions'],
            index=['Base Pairs Mapped to Intergenic Regions']
        ))
        region_plot = {'region_labels': ['Exonic Regions', 'Intronic Regions', 'Intergenic Regions'],
                'region_values': [Exonic_Regions, Intronic_Regions, Intergenic_Regions]}   
        self.add_content_item("data", STAR_plot=region_plot)

        self.stats.to_csv(self.stat_file, sep=':', header=False)

    @utils.add_log
    def ribo(self):
        human_ribo_fa = f'{ROOT_PATH}/data/rRNA/human_ribo.fasta'
        self.ribo_log = f'{self.outdir}/{self.sample}_ribo_log.txt'
        self.ribo_run_log = f'{self.outdir}/{self.sample}_ribo_run.log'
        cmd = (
            f'bbduk.sh '
            f'in1={self.fq} '
            f'ref={human_ribo_fa} '
            f'stats={self.ribo_log} '
            f'overwrite=t '
            f'> {self.ribo_run_log} 2>&1 '
        )
        Star.ribo.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
    
    @utils.add_log
    def STAR(self):
        cmd = [
            'STAR',
            '--runThreadN', str(self.thread),
            '--genomeDir', self.STAR_index,
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
        Star.STAR.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def picard(self):
        refFlat, _gtf, _ = glob_genomeDir(self.genomeDir)
        cmd = [
            'picard',
            '-Xmx20G',
            '-XX:ParallelGCThreads=4',
            'CollectRnaSeqMetrics',
            'I=%s' % (self.STAR_bam),
            'O=%s' % (self.picard_region_log),
            'REF_FLAT=%s' % (refFlat),
            'STRAND=NONE',
            'VALIDATION_STRINGENCY=SILENT']
        cmd_str = ' '.join(cmd)
        Star.picard.logger.info(cmd_str)
        subprocess.check_call(cmd)

    @utils.add_log
    def run(self):
        self.STAR()
        self.sort_bam()
        self.index_bam()
        self.picard()
        if self.debug:
            self.ribo()
        self.format_stat()
        self.clean_up()

    @utils.add_log
    def sort_bam(self):
        cmd = (
            f'samtools sort {self.unsort_STAR_bam} '
            f'-o {self.STAR_bam} '
            f'--threads {self.thread} '
        )
        Star.sort_bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def index_bam(self):
        cmd = f"samtools index {self.STAR_bam}"
        Star.index_bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)


def star(args):
    step_name = "star"
    runner = Star(args, step_name)
    runner.run()


def get_opts_star(parser, sub_program):
    parser.add_argument('--outFilterMatchNmin', help='STAR outFilterMatchNmin', default=0)
    parser.add_argument('--out_unmapped', help='out_unmapped', action='store_true')
    parser.add_argument('--genomeDir', help='genome directory')
    parser.add_argument('--STAR_param', help='STAR parameters', default="")
    parser.add_argument('--STAR_index', help='STAR index directory')
    parser.add_argument('--refFlat', help='refFlat file path')
    parser.add_argument('--outFilterMultimapNmax', help='STAR outFilterMultimapNmax', default=1)
    parser.add_argument('--starMem', help='starMem', default=30)
    if sub_program:
        parser.add_argument('--fq', required=True)
        parser.add_argument("--consensus_fq", action='store_true', help="input fastq is umi consensus")
        parser = s_common(parser)
