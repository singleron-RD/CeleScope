import argparse
import glob
import itertools
import os
from collections import defaultdict

import pandas as pd

import celescope
from celescope.tools.__init__ import FILTERED_MATRIX_DIR_SUFFIX
from celescope.tools import utils
from celescope.celescope import ArgFormatter
from celescope.__init__ import HELP_DICT

TOOLS_DIR = os.path.dirname(celescope.tools.__file__)
SAMPLE_TSV_REQUIRED_COLS = ['fastq_prefix', 'fastq_dir', 'sample_name']


def get_read(fastq_prefix, fastq_dirs, read='1'):
    read1_list = [f'_{read}', f'R{read}', f'R{read}_001']
    fq_list = ['fq', 'fastq']
    suffix_list = ["", ".gz"]
    fastq_dirs = fastq_dirs.split(",")
    read_pattern_list = [
        f'{fastq_dir}/{fastq_prefix}*{read}.{fq_str}{suffix}'
        for read in read1_list
        for fq_str in fq_list
        for suffix in suffix_list
        for fastq_dir in fastq_dirs
    ]
    fq_list = [glob.glob(read1_pattern) for read1_pattern in read_pattern_list]
    fq_list = (non_empty for non_empty in fq_list if non_empty)
    fq_list = sorted(list(itertools.chain(*fq_list)))
    if len(fq_list) == 0:
        print("Allowed R1 patterns:")
        for pattern in read_pattern_list:
            print(pattern)
        raise Exception(
            '\n'
            f'Invalid Read{read} path! \n'
            f'fastq_prefix: {fastq_prefix}\n'
            f'fastq_dir: {fastq_dirs}\n'
        )
    return fq_list


def get_fq(fastq_prefix, fastq_dir):
    """
    one (fastq_prefix, fastq_dir) combination can have multiple fq1 and fq2
    return fq1_list, fq2_list
    """
    fq1_list = get_read(fastq_prefix, fastq_dir, read='1')
    fq2_list = get_read(fastq_prefix, fastq_dir, read='2')
    if len(fq1_list) != len(fq2_list):
        raise Exception("Read1 and Read2 fastq number do not match!")

    return fq1_list, fq2_list


def get_fq_multiple(fastq_prefix_multiple, fastq_dir):
    fq1_list, fq2_list = [] , []
    fastq_prefix_list = fastq_prefix_multiple.split(',')
    for fastq_prefix in fastq_prefix_list:
        fastq_prefix = fastq_prefix.strip()
        fq1_list_tmp, fq2_list_tmp = get_fq(fastq_prefix, fastq_dir)
        fq1_list.extend(fq1_list_tmp)
        fq2_list.extend(fq2_list_tmp)

    fq1_str = ','.join(fq1_list)
    fq2_str = ','.join(fq2_list)
    
    return fq1_str, fq2_str

@utils.add_log
def parse_sample_tsv(sample_tsv):
    """
    Args:
        sample_tsv: sample tsv file with header.
    Returns:
        sample_dict:
        {'sample1':
            {'fastq_prefix': 'rna',
            'fastq_dir': '/celescope_test_data/rna/fastqs',
            'force_cell_num': 1000.0
        },
    """
    df = pd.read_csv(sample_tsv, sep='\t', header=0)
    df.fillna('None', inplace=True)

    for col in SAMPLE_TSV_REQUIRED_COLS:
        if not col in df.columns:
            raise Exception(f'{col} header is missing in mapfile! Required colnames: {col}')

    df = df.set_index('sample_name')
    sample_dict = df.to_dict('index')

    # add_fastq_path_to_dict
    for sample_name in sample_dict:
        fq1_str, fq2_str = get_fq_multiple(sample_dict[sample_name]['fastq_prefix'], sample_dict[sample_name]['fastq_dir'])
        sample_dict[sample_name]['fq1'] = fq1_str
        sample_dict[sample_name]['fq2'] = fq2_str

    return sample_dict


class Multi():

    def __init__(self, assay):
        self.__ASSAY__ = assay
        init_module = utils.find_assay_init(assay)
        self.__STEPS__ = init_module.__STEPS__
        self.__CONDA__ = os.path.basename(os.environ['CONDA_DEFAULT_ENV'])
        self.__APP__ = 'celescope'
        self.steps_not_run = ['mkref']

        # remove
        for step in self.steps_not_run:
            if step in self.__STEPS__:
                self.__STEPS__.remove(step)

        # add args
        self.parser = None
        self.common_args()
        self.step_args()

        # set
        self.args = None
        self.last_step = ''
        self.fq_suffix = ""
        self.steps_run = self.__STEPS__
        self.sample_dict = {}

        self.logdir = None
        self.sjm_dir = None
        self.sjm_file = None

        self.sjm_cmd = ''
        self.sjm_order = ''
        self.shell_dict = defaultdict(str)
  
        #snakemake
        self.snakemake_dir = None
        self.the_end_step = self.steps_run[-1]
        self.snakemake_shell_rule_all_dict = defaultdict(str)
        self.snakemake_shell_dict = defaultdict(str)

        self.outdir_dic = {}

    def common_args(self):
        readme = f'{self.__ASSAY__} multi-samples'
        parser = argparse.ArgumentParser(readme,
                                         formatter_class=ArgFormatter,
                                         conflict_handler='resolve')
        parser.add_argument(
            '--sample_tsv',
            help='''
Sample_tsv is a tab-delimited text file with header. Each line of mapfile represents paired-end fastq files. It contains at least 3 columns:
- `fastq_prefix`: Fastq file name prefix. Multiple fastq name prefix are separated by comma.
- `fastq_dir`: Fastq file directory path.
- `sample_name`: Unique name of the sample. It is the prefix of all output files corresponds to this sample.

Additional required columns:
- `matched_dir`: The single cell rna directory after running CeleScope is called `matched_dir`. 
This column in required in `tag`, `snp` and `capture_virus`, and optional in `vdj`.
- `background_snp`: Background snp file. It is required in `dynaseq`. 

Additional optional columns:
Any availabel arguments in `multi_{assay}` can be used in this file. 
- `expected_cell_num` Expected cell number.
- `force_cell_num`: Force celescope to use this number of cells. 
...

Example

Sample1 has 2 paired-end fastq files located in fastq_dir1. Sample2 has 1 paired-end fastq file located in fastq_dir2.
```
$cat ./sample.tsv
sample_name fastq_prefix    fastq_dir
sample1 fastq_prefix1	fastq_dir1,fastq_dir2
sample2 fastq_prefix3	fastq_dir3

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz

$ls fastq_dir2
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz

$ls fastq_dir2
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz
```
''',
            required=True)
        parser.add_argument('--mod', help='Which type of script to generate, `sjm` or `shell` or `snakemake`.',
            choices=['sjm', 'shell','snakemake'], default='sjm')
        parser.add_argument('--queue', help='Only works if the `--mod` selects `sjm`.')
        parser.add_argument('--rm_files', action='store_true',
            help='Remove redundant fastq and bam files after running.')
        parser.add_argument('--steps_run', 
            help='''
Steps to run. Multiple Steps are separated by comma. For example, if you only want to run `barcode` and `cutadapt`, 
use `--steps_run barcode,cutadapt`
''', 
            default='all')
        # sub_program parser do not have
        parser.add_argument('--outdir', help='Output directory.', default="./")
        parser.add_argument('--thread', help=HELP_DICT['thread'], default=4)
        parser.add_argument('--debug', help=HELP_DICT['debug'], action='store_true')
        self.parser = parser
        return parser

    def step_args(self):
        for step in self.__STEPS__:
            step_module = utils.find_step_module(self.__ASSAY__, step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            func_opts(self.parser, sub_program=False)



    def prepare(self):
        """
        parse_mapfile, make log dir, init script variables, init outdir_dic
        make sjm dir, sjm file
        """
        self.args = self.parser.parse_args()

        if self.args.gzip:
            self.fq_suffix = ".gz"
        if self.args.steps_run != 'all':
            self.steps_run = self.args.steps_run.strip().split(',')
        
        if self.args.mod == 'sjm':

            self.sjm_dir = f'{self.args.outdir}/sjm/'
            utils.check_mkdir(self.sjm_dir)
            self.logdir = self.args.outdir + '/log'
            utils.check_mkdir(self.logdir)

            self.sjm_file = f'{self.sjm_dir}/sjm.job'
            self.sjm_cmd = f'log_dir {self.logdir}\n'

        #add_snakemake
        if self.args.mod == 'snakemake':
            self.snakemake_dir = f'{self.args.outdir}/snakemake/'
            utils.check_mkdir(self.snakemake_dir)

        # parse_mapfile
        self.sample_dict = parse_sample_tsv(self.args.sample_tsv)
        

        for sample in self.sample_dict:
            self.outdir_dic[sample] = {}
            index = 0
            for step in self.__STEPS__:
                step_outdir = f"{self.args.outdir}/{sample}/{index:02d}.{step}"
                self.outdir_dic[sample].update({step: step_outdir})
                index += 1

    def generate_cmd(self, cmd, step, sample, m=1, x=1):
        if sample:
            sample = "_" + sample
        sched_options = f'sched_options -w n -cwd -V -l vf={m}g,p={x}'
        if self.args.queue:
            sched_options += f' -q {self.args.queue} '
        self.sjm_cmd += f'''
job_begin
    name {step}{sample}
    {sched_options}
    cmd source activate {self.__CONDA__}; {cmd}
job_end
'''


    def generate_snakemake_cmd(self, cmd, step, outdir, last_step, sample, x=1,fq=None):
        """
        this is snakemake need
        """
        if sample:
            sample = "_" + sample
        if step == 'sample' and fq != None:
            self.snakemake_shell_dict[sample] += f'''
rule {step}{sample}:
    input:
        {fq},
    output:
        directory("{outdir}"),
    threads:{x},
    params:
        step_name="celescope {self.__ASSAY__} {step}",
        cmd_line="{cmd}",
    message:"Executing {{params.step_name}} with {{threads}} threads.",
    conda:"{self.__CONDA__}",
    shell:
        "{{params.cmd_line}}"       
'''
        else:
            self.snakemake_shell_dict[sample] += f'''
rule {step}{sample}:
    input:
        rules.{last_step}{sample}.output,
    output:
        directory("{outdir}"),
    threads:{x},
    params:
        step_name="celescope {self.__ASSAY__} {step}",
        cmd_line="{cmd}",
    message:"Executing {{params.step_name}} with {{threads}} threads.",
    conda:"{self.__CONDA__}",
    shell:
        "{{params.cmd_line}}"       
'''


    def snakemake_rule_all(self,sample,outdir,step):
        """
        the rule all of snakemake
        """
        if sample:
            sample = "_" + sample
        if step == self.the_end_step:
            self.snakemake_shell_rule_all_dict[sample] += f'''
rule all:
    input:"{outdir}"      
'''         
        


    def process_cmd(self, cmd, step, sample, m=1, x=1):
        self.generate_cmd(cmd, step, sample, m=m, x=x)
        self.shell_dict[sample] += cmd + '\n'
        if self.last_step:
            self.sjm_order += f'order {step}_{sample} after {self.last_step}_{sample}\n'
        self.last_step = step



    def process_snakemake_cmd(self, cmd, step, outdir, sample, x=1,fq=None):
        """
        creat snakemake rule
        """
        self.generate_snakemake_cmd(cmd, step, outdir, self.last_step, sample, fq=fq, x=x) 
        self.snakemake_rule_all(sample,outdir,step)



    def get_args_from_multi_cli(self, step):
        """
        Get arguments from multi_{assay}
        """
        step_module = utils.find_step_module(self.__ASSAY__, step)
        func_opts = getattr(step_module, f"get_opts_{step}")
        step_parser = argparse.ArgumentParser(step_module)
        func_opts(step_parser, sub_program=False)
        known_args, unknown_args = step_parser.parse_known_args()
        return known_args, unknown_args

    def get_all_step_args(self, step):
        """
        Get all arguments from step cli
        """
        step_module = utils.find_step_module(self.__ASSAY__, step)
        func_opts = getattr(step_module, f"get_opts_{step}")
        step_parser = argparse.ArgumentParser(step_module)
        func_opts(step_parser, sub_program=True)
        all_step_args = [action.dest for action in step_parser._actions]
        if 'help' in all_step_args:
            all_step_args.remove('help')

        return all_step_args

    def get_cmd_line(self, step, sample):
        """ get cmd line 
            return str
        """
        known_args, _unknown_args = self.get_args_from_multi_cli(step)
        args_dict = known_args.__dict__
        step_prefix = (
            f'{self.__APP__} {self.__ASSAY__} {step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--thread {self.args.thread} '
        )
        cmd_line = step_prefix
        if self.args.debug:
            cmd_line += " --debug "

        # override all args value with value in sample_tsv
        all_step_args = self.get_all_step_args(step)
        for arg in all_step_args:
            if arg in self.sample_dict[sample]:
                value = self.sample_dict[sample][arg]
                if value != 'None':
                    cmd_line += f'--{arg} {value} '

        for arg in args_dict:
            # If in sample_tsv, override it
            if arg in self.sample_dict[sample]:
                print(f"WARNING: Command line argument `{arg}` of sample `{sample}` will be overrided by the value in sample_tsv")
            else:
                if args_dict[arg] is False:
                    continue
                if args_dict[arg] is True:
                    cmd_line += f'--{arg} '
                else:
                    if args_dict[arg]:
                        matches = [' ', '-']
                        arg_string = str(args_dict[arg])
                        if any(char in arg_string for char in matches):  # need quote
                            cmd_line += f'--{arg} "{arg_string}" '
                        else:
                            cmd_line += f'--{arg} {arg_string} '

        return cmd_line

    def sample(self, sample):
        step = "sample"
        cmd_line = self.get_cmd_line(step, sample)
        fq = self.sample_dict[sample]['fq1'].split(",")
        out_dir = f'{self.outdir_dic[sample][step]}'
        cmd = (
            f'{cmd_line} '
        )

        self.process_snakemake_cmd(cmd, step, out_dir,sample,fq=fq,x=1)  
        self.process_cmd(cmd, step, sample, m=1, x=1)
       

    def barcode(self, sample):
        step = "barcode"
        cmd_line = self.get_cmd_line(step, sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        cmd = (
            f'{cmd_line} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=5, x=1)  
        

    def cutadapt(self, sample):
        step = "cutadapt"
        fq = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq{self.fq_suffix}'
        out_dir = f'{self.outdir_dic[sample][step]}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=5, x=1)
       

    def star(self, sample):
        step = 'star'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        out_dir = f'{self.outdir_dic[sample][step]}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=self.args.thread)
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread) 
        

    def featureCounts(self, sample):
        step = 'featureCounts'
        input_bam = f'{self.outdir_dic[sample]["star"]}/{sample}_Aligned.sortedByCoord.out.bam'
        out_dir = f'{self.outdir_dic[sample][step]}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--input {input_bam} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=self.args.thread)
        self.process_cmd(cmd, step, sample, m=5, x=self.args.thread)
        

    def count(self, sample):
        step = 'count'
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_name_sorted.bam'
        out_dir = f'{self.outdir_dic[sample][step]}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=10, x=1)
        

    def analysis(self, sample):
        step = 'analysis'
        matrix_file = f'{self.outdir_dic[sample]["count"]}/{sample}_{FILTERED_MATRIX_DIR_SUFFIX[0]}'
        out_dir = f'{self.outdir_dic[sample][step]}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--matrix_file {matrix_file} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=10, x=1)
        

    def consensus(self, sample):
        step = 'consensus'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        out_dir = f'{self.outdir_dic[sample][step]}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=5, x=1)    
        outfile = f'{self.outdir_dic[sample][step]}/{sample}_consensus.fq'
        return outfile

    def run_steps(self):
        for sample in self.sample_dict:
            self.last_step = ''
            for step in self.steps_run:
                if step in self.steps_not_run:
                    continue
                try:
                    method_to_call = getattr(self, step)
                except AttributeError as attr_not_exist:
                    raise NotImplementedError(
                        "Class `{}` does not implement `{}`".format(self.__class__.__name__, step)
                    ) from attr_not_exist
                method_to_call(sample)

    def merge_report(self):
        step = "merge_report"
        steps_str = ",".join(self.__STEPS__)
        samples = ','.join(self.sample_dict.keys())
        app = TOOLS_DIR + '/merge_table.py'
        cmd = (
            f'python {app} --samples {samples} '
            f'--steps {steps_str} --outdir {self.args.outdir}'
        )
        if self.args.rm_files:
            cmd += ' --rm_files'
        self.generate_cmd(cmd, step, sample="")
        for sample in self.sample_dict:
            self.sjm_order += f'order {step} after {self.last_step}_{sample}\n'

    def end(self):
        if self.args.mod == 'sjm':
            self.merge_report()
            with open(self.sjm_file, 'w') as fh:
                fh.write(self.sjm_cmd + '\n')
                fh.write(self.sjm_order)
        if self.args.mod == 'shell':
            os.system('mkdir -p ./shell/')
            for sample in self.shell_dict:
                with open(f'./shell/{sample}.sh', 'w') as f:
                    f.write(self.shell_dict[sample])
        if self.args.mod == "snakemake":
            for sample in self.snakemake_shell_dict:
                with open(f'{self.snakemake_dir}/Snakefile{sample}','w') as fh:
                    fh.write(self.snakemake_shell_rule_all_dict[sample] + '\n')
                    fh.write(self.snakemake_shell_dict[sample])

    def run(self):
        self.prepare()
        self.run_steps()
        self.end()


