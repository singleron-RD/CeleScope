import os
import glob
import sys
import argparse
import re
from collections import defaultdict
from celescope.tools.utils import *
from celescope.tools.__init__ import __PATTERN_DICT__

class Multi():

    def __init__(self, assay):
        self.__ASSAY__ = assay
        init_module = find_assay_init(assay)
        self.__STEPS__ = init_module.__STEPS__
        self.__CONDA__ = os.environ['CONDA_DEFAULT_ENV']
        self.__APP__ = 'celescope'
        self.col4_default = None
        self.last_step = ''
        self.step_prefix = (
            '{self.__APP__} {self.__ASSAY__} {step}'
            '--outdir {self.outdir_dic[sample][step]} '
            '--sample {sample} '
            '--assay {self.__ASSAY__} '
        )
        self.args = None

    def common_args(self):
        readme = f'{self.__ASSAY__} multi-samples'
        parser = argparse.ArgumentParser(readme, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
        parser.add_argument('--rm_files', action='store_true', help='remove redundant fq.gz and bam after running')
        parser.add_argument('--steps_run', help='steps to run', default='all')
        # sub_program parser do not have
        parser.add_argument('--outdir', help='output dir', default="./")
        parser.add_argument('--debug', help='debug or not', action='store_true')
        parser.add_argument('--thread', help='thread', default=4)
        self.parser = parser
        return parser

    def read_common_args(self):
        self.not_gzip_str = Multi.arg_str(self.args.not_gzip, 'not_gzip')
        if os.path.basename(self.__CONDA__) == 'celescope_RD':
            self.debug_str = '--debug'
        else:
            self.debug_str = Multi.arg_str(self.args.debug, 'debug')

    def step_args(self):
        for step in self.

    def parse_step_args(self, step):
        step_module = find_step_module(self.__ASSAY__, step)
        func_opts = getattr(step_module, f"get_opts_{step}")
        parser = func_opts(self.parser, sub_program=False)
        args = parser.parse_known_args()
        return args


    def parse_args(self):
        self.common_args()
        self.step_args()
        self.args = self.parser.parse_args()
        if not self.args.not_gzip:
            self.fq_suffix = ".gz"
        else:
            self.fq_suffix = ""


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
            f"{f'{self.step_prefix}'} "
            f'--chemistry {self.args.chemistry} '
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
            f'--chemistry {self.args.chemistry} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} '
            f'--pattern {self.args.pattern} --whitelist {self.args.whitelist} --linker {self.args.linker} '
            f'--lowQual {self.args.lowQual} --thread {self.args.thread} '
            f'--lowNum {self.args.lowNum} '
            f'{self.args.allowNoPolyT_str} '
            f'{self.allowNoLinker_str} '
            f'{self.noLinker_str} '
            f'{self.nopolyT_str} '
            f'{self.not_gzip_str} '
            f'--probe_file {self.probe_file} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def cutadapt(self, sample):
        # adapt
        step = "cutadapt"
        fq = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq{self.fq_suffix}'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step } '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--overlap {self.overlap} '
            f'--minimum_length {self.minimum_length} '
            f'--insert {self.insert} '
            f'--adapter_fasta {self.adapter_fasta} '
            f'{self.not_gzip_str} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def STAR(self, sample):
        step = 'STAR'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
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
            f'--STAR_param \"{self.STAR_param}\" '
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
            f'--force_cell_num {self.col4_dict[sample]} '
            f'--genomeDir {self.genomeDir} '
            f'--cell_calling_method {self.cell_calling_method} '
            f'--expected_cell_num {self.expected_cell_num} '
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

    def consensus_args(self):
        self.parser.add_argument("--threshold", help='valid base threshold', default=0.5)

    def read_consensus_args(self):
        self.threshold = self.args.threshold

    def consensus(self, sample):
        step = 'consensus'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step } '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--threshold {self.threshold} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)
        outfile = f'{self.outdir_dic[sample][step]}/{sample}_consensus.fq'
        return outfile

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
        self.parse_args()
        self.prepare()
        self.run_steps()
        self.end()



    