import os
import glob
import sys
import argparse
import re
from collections import defaultdict
from celescope.__init__ import __CONDA__
from celescope.tools.utils import merge_report
from celescope.tools.utils import parse_map_col4, link_data
from celescope.tools.__init__ import __PATTERN_DICT__

class Multi():

    def __init__(self, __ASSAY__, __STEPS__, __CONDA__):
        self.__ASSAY__ = __ASSAY__
        self.__STEPS__ = __STEPS__
        self.__CONDA__ = __CONDA__
        self.__APP__ = 'celescope'
        self.col4_default = 'auto'
        self.last_step = ''

    def multi_opts(self):
        readme = f'{self.__ASSAY__} multi-samples'
        parser = argparse.ArgumentParser(readme)
        parser.add_argument('--mod', help='mod, sjm or shell', choices=['sjm', 'shell'], default='sjm')
        parser.add_argument(
            '--mapfile',
            help='''
                tsv file, 4 columns:
                1st col: LibName;
                2nd col: DataDir;
                3rd col: SampleName;
                4th col: Cell number or match_dir, optional;
            ''',
            required=True)
        parser.add_argument('--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version', default='auto')
        parser.add_argument('--whitelist', help='cellbarcode list')
        parser.add_argument('--linker', help='linker')
        parser.add_argument('--pattern', help='read1 pattern')
        parser.add_argument('--outdir', help='output dir', default="./")
        parser.add_argument(
            '--adapt',
            action='append',
            help='adapter sequence',
            default=[
                'polyT=A{15}',
                'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'])
        parser.add_argument(
            '--minimum_length',
            dest='minimum_length',
            help='minimum_length',
            default=20)
        parser.add_argument(
            '--nextseq-trim',
            dest='nextseq_trim',
            help='nextseq_trim',
            default=20)
        parser.add_argument('--overlap', help='minimum overlap length', default=10)
        parser.add_argument('--lowQual', type=int, help='max phred of base as lowQual', default=0)
        parser.add_argument(
            '--lowNum',
            type=int,
            help='max number with lowQual allowed',
            default=2)
        parser.add_argument(
            '--rm_files',
            action='store_true',
            help='remove redundant fq.gz and bam after running')
        parser.add_argument('--steps_run', help='steps to run', default='all')
        parser.add_argument('--debug', help='debug or not', action='store_true')
        self.parser = parser
        return parser

    def STAR_args(self):
        self.parser.add_argument('--starMem', help='starMem', default=30)
        self.parser.add_argument('--genomeDir', help='genome index dir', required=True)
        self.parser.add_argument(
            '--gtf_type',
            help='Specify attribute type in GTF annotation, default=exon',
            default='exon')
        self.parser.add_argument('--thread', help='thread', default=6)
        self.parser.add_argument('--out_unmapped', help='out_unmapped', action='store_true')
        self.parser.add_argument('--outFilterMatchNmin', help='STAR outFilterMatchNmin', default=0)

    def analysis_args(self):
        self.parser.add_argument('--save_rds', action='store_true', help='write rds to disk')
        self.parser.add_argument('--type_marker_tsv', help='cell type marker tsv')

    def custome_args(self):
        self.STAR_args()
        self.analysis_args()

    def parse_args(self):
        self.multi_opts()
        self.custome_args()
        self.args = self.parser.parse_args()
        # read args
        self.outdir = self.args.outdir
        self.chemistry = self.args.chemistry
        self.pattern = self.args.pattern
        self.whitelist = self.args.whitelist
        self.linker = self.args.linker
        self.lowQual = self.args.lowQual
        self.lowNum = self.args.lowNum
        self.overlap = self.args.overlap
        self.mod = self.args.mod
        self.rm_files = self.args.rm_files
        self.steps_run = self.args.steps_run
        if self.__CONDA__ == 'celescope_RD':
            self.debug_str = '--debug'
        else:
            self.debug_str = Multi.arg_str(self.args.debug, 'debug')

    @staticmethod
    def arg_str(arg, arg_name):
        '''
        return action store_true arguments as string
        '''
        if arg:
            return '--' + arg_name
        return ''

    def read_STAR_args(self):
        self.thread = self.args.thread
        self.genomeDir = self.args.genomeDir
        self.starMem = self.args.starMem
        self.gtf_type = self.args.gtf_type
        self.out_unmapped = Multi.arg_str(self.args.out_unmapped, 'out_unmapped')
        self.outFilterMatchNmin = self.args.outFilterMatchNmin
    
    def read_analysis_args(self):
        self.save_rds = self.args.save_rds
        self.save_rds_str = Multi.arg_str(self.save_rds, 'save_rds')
        self.type_marker_tsv = self.args.type_marker_tsv

    def read_custome_args(self):
        self.read_STAR_args()
        self.read_analysis_args()

    def prepare(self):
        # parse_mapfile
        self.fq_dict, self.col4_dict = parse_map_col4(self.args.mapfile, self.col4_default)

        # link
        link_data(self.outdir, self.fq_dict)      

        # mk log dir
        self.logdir = self.outdir + '/log'
        os.system('mkdir -p %s' % (self.logdir))

        # script init
        self.sjm_cmd = 'log_dir %s\n' % (self.logdir)
        self.sjm_order = ''
        self.shell_dict = defaultdict(str)

        # outdir dict
        self.outdir_dic = {}
        for sample in self.fq_dict:
            self.outdir_dic[sample] = {}
            index = 0
            for step in self.__STEPS__:
                step_outdir = f"{self.outdir}/{sample}/{index:02d}.{step}"
                self.outdir_dic[sample].update({step: step_outdir})
                index += 1
    
    def generate_cmd(self, cmd, step, sample, m, x):
        self.sjm_cmd += f'''
job_begin
    name {step}_{sample}
    sched_options -w n -cwd -V -l vf={m}g,p={x}
    cmd source activate {self.__CONDA__}; {cmd}
job_end
'''

    def process_cmd(self, cmd, step, sample, m=1, x=1):
        self.generate_cmd(cmd, step, sample, m=m, x=x)
        self.shell_dict[sample] += cmd + '\n'
        if self.last_step:
            self.sjm_order += f'order {step}_{sample} after {self.last_step}_{sample}\n'
        self.last_step = step

    def generate_first(self, cmd, step, sample, m=1, x=1):
        self.generate_cmd(cmd, step, sample, m=1, x=1)
        self.shell_dict[sample] += cmd + '\n'
        self.last_step = step

    def generate_other(self, cmd, step, sample, m=1, x=1):
        self.generate_cmd(cmd, step, sample, m=1, x=1)
        self.shell_dict[sample] += cmd + '\n'
        self.sjm_order += f'order {step}_{sample} after {self.last_step}_{sample}\n'
        self.last_step = step

    def sample(self, sample):
        step = "sample"
        arr = self.fq_dict[sample]
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--chemistry {self.chemistry} '
            f'--fq1 {arr[0]}'
        )
        self.process_cmd(cmd, step, sample, m=1, x=1)
    
    def barcode(self, sample):
        # barcode
        arr = self.fq_dict[sample]
        step = "barcode"
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--chemistry {self.chemistry} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} '
            f'--pattern {self.pattern} --whitelist {self.whitelist} --linker {self.linker} '
            f'--lowQual {self.lowQual} --thread {self.thread} '
            f'--lowNum {self.lowNum} '

        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def cutadapt(self, sample):
        # adapt
        step = "cutadapt"
        fq = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq.gz'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step } '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--overlap {self.overlap} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def STAR(self, sample):
        step = 'STAR'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--genomeDir {self.genomeDir} '
            f'--thread {self.thread} '
            f'{self.debug_str} '
            f'--outFilterMatchNmin {self.outFilterMatchNmin} '
            f'{self.out_unmapped} '
        )
        self.process_cmd(cmd, step, sample, m=self.starMem, x=self.thread)

    def featureCounts(self, sample):
        step = 'featureCounts'
        input = f'{self.outdir_dic[sample]["STAR"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--input {input} --gtf_type {self.gtf_type} '
            f'--genomeDir {self.genomeDir} '
            f'--thread {self.thread} '
            f'--gtf_type {self.gtf_type} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=self.thread)

    def count(self, sample):
        step = 'count'
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_name_sorted.bam'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--bam {bam} '
            f'--cells {self.col4_dict[sample]} '
            f'--genomeDir {self.genomeDir} '
        )
        self.process_cmd(cmd, step, sample, m=10, x=1)

    def analysis(self, sample):
        step = 'analysis'
        matrix_file = f'{self.outdir_dic[sample]["count"]}/{sample}_matrix.tsv.gz'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step } '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--matrix_file {matrix_file} '
            f'{self.save_rds_str} '
            f'--type_marker_tsv {self.type_marker_tsv} '
        )
        self.process_cmd(cmd, step, sample, m=10, x=1)

    def run_steps(self):
        if self.steps_run == 'all':
            self.steps_run = self.__STEPS__
        elif self.steps_run:
            self.steps_run = self.steps_run.strip().split(',')

        for sample in self.fq_dict:
            self.last_step = ''
            for step in self.steps_run:
                eval(f'self.{step}(sample)')

    def end(self):
        if self.mod == 'sjm':
            step = 'merge_report'
            merge_report(
                self.fq_dict,
                self.__STEPS__,
                self.last_step,
                self.sjm_cmd,
                self.sjm_order,
                self.logdir,
                self.__CONDA__,
                self.outdir,
                self.rm_files,
            )
        if self.mod == 'shell':
            os.system('mkdir -p ./shell/')
            for sample in self.shell_dict:
                with open(f'./shell/{sample}.sh', 'w') as f:
                    f.write(self.shell_dict[sample])

    def run(self):
        self.multi_opts()
        self.custome_args()
        self.parse_args()
        self.read_custome_args()
        self.prepare()
        self.run_steps()
        self.end()



    