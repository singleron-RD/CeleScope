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
        self.__CONDA__ = os.path.basename(os.environ['CONDA_DEFAULT_ENV'])
        self.__APP__ = 'celescope'
        self.col4_default = None
        self.last_step = ''
        self.args = None

    def common_args(self):
        readme = f'{self.__ASSAY__} multi-samples'
        parser = argparse.ArgumentParser(readme, 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            conflict_handler='resolve')
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

    def step_args(self):
        for step in self.__STEPS__:
            step_module = find_step_module(self.__ASSAY__, step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            func_opts(self.parser, sub_program=False)

    def parse_args(self):
        self.common_args()
        self.step_args()
        self.args = self.parser.parse_args()
        if not self.args.not_gzip:
            self.fq_suffix = ".gz"
        else:
            self.fq_suffix = ""

    @staticmethod
    def parse_map_col4(mapfile, default_val):
        fq_dict = defaultdict(list)
        col4_dict = {}
        with open(mapfile) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                tmp = line.split()
                library_id, library_path, sample_name = tmp[:3]
                if len(tmp) == 4:
                    col4 = tmp[3]
                else:
                    col4 = default_val
                fq1, fq2 = get_fq(library_id, library_path)

                if sample_name in fq_dict:
                    fq_dict[sample_name][0].append(fq1)
                    fq_dict[sample_name][1].append(fq2)
                else:
                    fq_dict[sample_name] = [[fq1], [fq2]]
                    col4_dict[sample_name] = col4


        for sample_name in fq_dict:
            fq_dict[sample_name][0] = ",".join(fq_dict[sample_name][0])
            fq_dict[sample_name][1] = ",".join(fq_dict[sample_name][1])

        if not fq_dict:
            raise Exception('empty mapfile!')
        return fq_dict, col4_dict


    def prepare(self):
        # parse_mapfile
        self.fq_dict, self.col4_dict = self.parse_map_col4(self.args.mapfile, self.col4_default)

        # link
        link_data(self.args.outdir, self.fq_dict)      

        # mk log dir
        self.logdir = self.args.outdir + '/log'
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
                step_outdir = f"{self.args.outdir}/{sample}/{index:02d}.{step}"
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

    def parse_step_args(self, step):
        step_module = find_step_module(self.__ASSAY__, step)
        func_opts = getattr(step_module, f"get_opts_{step}")
        step_parser = argparse.ArgumentParser(step_module)
        func_opts(step_parser, sub_program=False)
        args = step_parser.parse_known_args()
        return args

    def get_cmd_line(self, step, sample):
        """ get cmd line without input
            return str
        """
        args = self.parse_step_args(step)
        args_dict = args[0].__dict__
        step_prefix = (
            f'{self.__APP__} {self.__ASSAY__} {step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--thread {self.args.thread} '
        )
        cmd_line = step_prefix
        if self.__CONDA__ == "celescope_RD":
            cmd_line += " --debug "
        for arg in args_dict:
            if args_dict[arg] is False:
                continue
            if args_dict[arg] is True:
                cmd_line += f'--{arg} '
            else:
                if args_dict[arg]:
                    matches = [' ', '-']
                    arg_string = str(args_dict[arg])
                    if any(char in arg_string for char in matches): # need quote
                        cmd_line += f'--{arg} "{arg_string}" '
                    else:
                        cmd_line += f'--{arg} {arg_string} '

        return cmd_line

    def sample(self, sample):
        step = "sample"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr[0]} '
        )
        self.process_cmd(cmd, step, sample, m=1, x=1)
    
    def barcode(self, sample):
        step = "barcode"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def cutadapt(self, sample):
        step = "cutadapt"
        fq = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq{self.fq_suffix}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def STAR(self, sample):
        step = 'STAR'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def featureCounts(self, sample):
        step = 'featureCounts'
        input = f'{self.outdir_dic[sample]["STAR"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--input {input} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=self.args.thread)

    def count(self, sample):
        step = 'count'
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_name_sorted.bam'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--force_cell_num {self.col4_dict[sample]} '
        )

        self.process_cmd(cmd, step, sample, m=10, x=1)

    def analysis(self, sample):
        step = 'analysis'
        matrix_file = f'{self.outdir_dic[sample]["count"]}/{sample}_matrix.tsv.gz'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--matrix_file {matrix_file} '
        )
        self.process_cmd(cmd, step, sample, m=10, x=1)

    def consensus(self, sample):
        step = 'consensus'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)
        outfile = f'{self.outdir_dic[sample][step]}/{sample}_consensus.fq'
        return outfile

    def run_steps(self):
        if self.args.steps_run == 'all':
            self.steps_run = self.__STEPS__
        elif self.args.steps_run:
            self.steps_run = self.args.steps_run.strip().split(',')

        for sample in self.fq_dict:
            self.last_step = ''
            for step in self.steps_run:
                eval(f'self.{step}(sample)')

    def end(self):
        if self.args.mod == 'sjm':
            step = 'merge_report'
            merge_report(
                self.fq_dict,
                self.__STEPS__,
                self.last_step,
                self.sjm_cmd,
                self.sjm_order,
                self.logdir,
                self.__CONDA__,
                self.args.outdir,
                self.args.rm_files,
            )
        if self.args.mod == 'shell':
            os.system('mkdir -p ./shell/')
            for sample in self.shell_dict:
                with open(f'./shell/{sample}.sh', 'w') as f:
                    f.write(self.shell_dict[sample])

    def run(self):
        self.parse_args()
        self.prepare()
        self.run_steps()
        self.end()



    