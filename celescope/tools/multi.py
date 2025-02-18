import argparse
import glob
import itertools
import os
import sys
from collections import defaultdict

import celescope
from celescope.tools.__init__ import FILTERED_MATRIX_DIR_SUFFIX, STAR_BAM_SUFFIX
from celescope.tools import utils
from celescope.celescope import ArgFormatter
from celescope.__init__ import HELP_DICT
from celescope.tools.make_ref import MakeRef

TOOLS_DIR = os.path.dirname(celescope.tools.__file__)


class Multi:
    def __init__(self, assay):
        self.__ASSAY__ = assay
        init_module = utils.find_assay_init(assay)
        self.STEPS = init_module.STEPS
        self.REMOVE_FROM_MULTI = getattr(init_module, "REMOVE_FROM_MULTI", set())
        self.REMOVE_FROM_MULTI.add("mkref")
        try:
            self.__CONDA__ = os.path.basename(os.environ["CONDA_DEFAULT_ENV"])
        except KeyError:
            print("CONDA_DEFAULT_ENV is not set. sjm mode may not available.")
            self.__CONDA__ = "celescope"
        self.__APP__ = "celescope"

        # remove
        for step in self.REMOVE_FROM_MULTI:
            if step in self.STEPS:
                self.STEPS.remove(step)

        # add args
        self.parser = None
        self.common_args()
        self.step_args()

        # set
        self.args = None
        self.col4_default = None
        self.last_step = ""
        self.steps_run = self.STEPS
        self.fq_dict = {}
        self.col4_dict = {}
        self.col5_dict = {}
        self.logdir = None
        self.sjm_dir = None
        self.sjm_file = None

        self.sjm_cmd = ""
        self.sjm_order = ""
        self.shell_dict = defaultdict(str)

        self.outdir_dic = {}

    def common_args(self):
        readme = f"{self.__ASSAY__} multi-samples"
        parser = argparse.ArgumentParser(
            readme, formatter_class=ArgFormatter, conflict_handler="resolve"
        )
        parser.add_argument(
            "--mapfile",
            help="""
Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  
4th column: The 4th column has different meaning for each assay. The single cell rna directory after running CeleScope is called `matched_dir`.

- `vdj` Required, matched_dir.
- `tag` Required, matched_dir.
- `dynaseq` Optional, forced cell number.
- `snp` Required, matched_dir.
- `capture_virus` Required, matched_dir.
- `fusion` Required, matched_dir.
- `citeseq` Required, matched_dir.
- `flv_trust4` Required, matched_dir.
- `sweetseq` Required, matched_dir.
 
5th column:
- `dynaseq` Required, background snp file.

Example

Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.
```
$cat ./my.mapfile
fastq_prefix1	fastq_dir1	sample1
fastq_prefix2	fastq_dir2	sample1
fastq_prefix3	fastq_dir1	sample2

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
```
""",
            required=True,
        )
        parser.add_argument(
            "--mod",
            help="Which type of script to generate, `sjm` or `shell`.",
            choices=["sjm", "shell"],
            default="sjm",
        )
        parser.add_argument("--queue", help="Only works if the `--mod` selects `sjm`.")
        parser.add_argument(
            "--rm_files",
            action="store_true",
            help="Remove redundant fastq and bam files after running.",
        )
        parser.add_argument(
            "--steps_run",
            help="""
Steps to run. Multiple Steps are separated by comma. For example, if you only want to run `barcode` and `cutadapt`, 
use `--steps_run barcode,cutadapt`
""",
            default="all",
        )
        # sub_program parser do not have
        parser.add_argument("--outdir", help="Output directory.", default="./")
        parser.add_argument("--thread", help=HELP_DICT["thread"], default=16)
        parser.add_argument("--debug", help=HELP_DICT["debug"], action="store_true")
        self.parser = parser
        return parser

    def step_args(self):
        for step in self.STEPS:
            step_module = utils.find_step_module(self.__ASSAY__, step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            func_opts(self.parser, sub_program=False)

    @staticmethod
    @utils.add_log
    def parse_mapfile(mapfile, default_val):
        fq_dict = defaultdict(lambda: defaultdict(list))
        col4_dict = {}
        col5_dict = {}
        with open(mapfile) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                line_split = line.split()
                library_id, library_path, sample_name = line_split[:3]
                if len(line_split) >= 4:
                    col4 = line_split[3]
                else:
                    col4 = default_val
                fq1, fq2 = get_fq(library_id, library_path)

                fq_dict[sample_name]["fq1"].append(fq1)
                fq_dict[sample_name]["fq2"].append(fq2)
                fq_dict[sample_name]["col4"].append(col4)
                col4_dict[sample_name] = col4
                if len(line_split) == 5:
                    col5_dict[sample_name] = line_split[4]

        for sample_name in fq_dict:
            fq_dict[sample_name]["fq1_str"] = ",".join(fq_dict[sample_name]["fq1"])
            fq_dict[sample_name]["fq2_str"] = ",".join(fq_dict[sample_name]["fq2"])

        if not fq_dict:
            raise Exception("empty mapfile!")
        return fq_dict, col4_dict, col5_dict

    def link_data(self):
        raw_dir = f"{self.args.outdir}/data_give/rawdata"
        os.system("mkdir -p %s" % (raw_dir))
        with open(raw_dir + "/ln.sh", "w") as fh:
            fh.write("cd %s\n" % (raw_dir))
            for s, arr in self.fq_dict.items():
                fh.write("ln -sf %s %s\n" % (arr["fq1_str"], s + "_1.fq.gz"))
                fh.write("ln -sf %s %s\n" % (arr["fq2_str"], s + "_2.fq.gz"))

    def prepare(self):
        """
        parse_mapfile, make log dir, init script variables, init outdir_dic
        make sjm dir, sjm file
        """
        self.args = self.parser.parse_args()
        self.check_genome()

        if self.args.steps_run != "all":
            self.steps_run = self.args.steps_run.strip().split(",")

        if self.args.mod == "sjm":
            self.sjm_dir = f"{self.args.outdir}/sjm/"
            utils.check_mkdir(self.sjm_dir)
            self.logdir = self.args.outdir + "/log"
            utils.check_mkdir(self.logdir)

            self.sjm_file = f"{self.sjm_dir}/sjm.job"
            self.sjm_cmd = f"log_dir {self.logdir}\n"

        # parse_mapfile
        self.fq_dict, self.col4_dict, self.col5_dict = self.parse_mapfile(
            self.args.mapfile, self.col4_default
        )

        for sample in self.fq_dict:
            self.outdir_dic[sample] = {}
            index = 0
            for step in self.STEPS:
                step_outdir = f"{self.args.outdir}/{sample}/{index:02d}.{step}"
                self.outdir_dic[sample].update({step: step_outdir})
                index += 1
            # add outs dir
            self.outdir_dic[sample]["outs"] = f"{self.args.outdir}/{sample}/outs"

    def generate_cmd(self, cmd, step, sample, m=1, x=1):
        if sample:
            sample = "_" + sample
        sched_options = f"sched_options -w n -cwd -V -l vf={m}g,p={x}"
        if self.args.queue:
            sched_options += f" -q {self.args.queue} "
        self.sjm_cmd += f"""
job_begin
    name {step}{sample}
    {sched_options}
    cmd source activate {self.__CONDA__}; {cmd}
job_end
"""

    def process_cmd(self, cmd, step, sample, m=1, x=1):
        self.generate_cmd(cmd, step, sample, m=m, x=x)
        self.shell_dict[sample] += cmd + "\n"
        if self.last_step:
            self.sjm_order += f"order {step}_{sample} after {self.last_step}_{sample}\n"
        self.last_step = step

    def parse_step_args(self, step):
        step_module = utils.find_step_module(self.__ASSAY__, step)
        func_opts = getattr(step_module, f"get_opts_{step}")
        step_parser = argparse.ArgumentParser()
        func_opts(step_parser, sub_program=False)
        args = step_parser.parse_known_args()
        return args

    def get_cmd_line(self, step, sample):
        """get cmd line without input
        return str
        """
        args = self.parse_step_args(step)
        args_dict = args[0].__dict__
        step_prefix = (
            f"{self.__APP__} {self.__ASSAY__} {step} "
            f"--outdir {self.outdir_dic[sample][step]} "
            f"--sample {sample} "
            f"--thread {self.args.thread} "
        )
        cmd_line = step_prefix
        if self.args.debug:
            cmd_line += " --debug "
        for arg in args_dict:
            if args_dict[arg] is False:
                continue
            if args_dict[arg] is True:
                cmd_line += f"--{arg} "
            else:
                if args_dict[arg]:
                    matches = [" ", "-"]
                    arg_string = str(args_dict[arg])
                    if any(char in arg_string for char in matches):  # need quote
                        cmd_line += f'--{arg} "{arg_string}" '
                    else:
                        cmd_line += f"--{arg} {arg_string} "

        return cmd_line

    def sample(self, sample):
        step = "sample"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} '
        self.process_cmd(cmd, step, sample, m=1, x=1)

    def barcode(self, sample):
        step = "barcode"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def cutadapt(self, sample):
        step = "cutadapt"
        fq = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--fq {fq} "
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def star(self, sample):
        step = "star"
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--fq {fq} "
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def prep_map(self, sample):
        step = "prep_map"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def featureCounts(self, sample):
        step = "featureCounts"
        if "star" in self.STEPS:
            prev = "star"
        elif "prep_map" in self.STEPS:
            prev = "prep_map"
        else:
            sys.exit("To use featureCounts, star or prep must in the steps!")
        input_bam = f"{self.outdir_dic[sample][prev]}/{sample}_{STAR_BAM_SUFFIX}"
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--input {input_bam} "
        self.process_cmd(cmd, step, sample, m=5, x=self.args.thread)

    def count(self, sample):
        step = "count"
        count_detail = (
            f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_count_detail.txt'
        )
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f"{cmd_line} "
            f"--count_detail {count_detail} "
            f"--force_cell_num {self.col4_dict[sample]} "
        )

        self.process_cmd(cmd, step, sample, m=20, x=1)

    def analysis(self, sample):
        step = "analysis"
        matrix_file = (
            f'{self.outdir_dic[sample]["count"]}/{sample}_{FILTERED_MATRIX_DIR_SUFFIX}'
        )
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--matrix_file {matrix_file} "
        self.process_cmd(cmd, step, sample, m=10, x=1)

    def consensus(self, sample):
        step = "consensus"
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--fq {fq} "
        self.process_cmd(cmd, step, sample, m=5, x=1)
        outfile = f"{self.outdir_dic[sample][step]}/{sample}_consensus.fq"
        return outfile

    def run_steps(self):
        for sample in self.fq_dict:
            self.last_step = ""
            for step in self.steps_run:
                try:
                    method_to_call = getattr(self, step)
                except AttributeError as attr_not_exist:
                    raise NotImplementedError(
                        "Class `{}` does not implement `{}`".format(
                            self.__class__.__name__, step
                        )
                    ) from attr_not_exist
                method_to_call(sample)

    def merge_report(self):
        step = "merge_report"
        steps_str = ",".join(self.STEPS)
        steps_str = steps_str.replace("prep_map", "barcode,cutadapt,star")
        steps_str = steps_str.replace("starsolo", "demultiplexing,mapping,cells")
        samples = ",".join(self.fq_dict.keys())
        app = TOOLS_DIR + "/merge_table.py"
        cmd = (
            f"python {app} --samples {samples} "
            f"--steps {steps_str} --outdir {self.args.outdir}"
        )
        if self.args.rm_files:
            cmd += " --rm_files"
        self.generate_cmd(cmd, step, sample="")
        for sample in self.fq_dict:
            self.sjm_order += f"order {step} after {self.last_step}_{sample}\n"

    def end(self):
        if self.args.mod == "sjm":
            self.merge_report()
            with open(self.sjm_file, "w") as fh:
                fh.write(self.sjm_cmd + "\n")
                fh.write(self.sjm_order)
        if self.args.mod == "shell":
            os.system("mkdir -p ./shell/")
            for sample in self.shell_dict:
                with open(f"./shell/{sample}.sh", "w") as f:
                    f.write(self.shell_dict[sample])

    def check_genome(self):
        for arg, val in vars(self.args).items():
            if val and "genomeDir" in arg:
                MakeRef.get_config(val)

    def run(self):
        self.prepare()
        self.run_steps()
        self.end()


def get_read(library_id, library_path, read="1"):
    read1_list = [f"_{read}", f"R{read}", f"R{read}_001"]
    fq_list = ["fq", "fastq"]
    suffix_list = ["", ".gz"]
    read_pattern_list = [
        f"{library_path}/{library_id}*{read}.{fq_str}{suffix}"
        for read in read1_list
        for fq_str in fq_list
        for suffix in suffix_list
    ]
    fq_list = [glob.glob(read1_pattern) for read1_pattern in read_pattern_list]
    fq_list = (non_empty for non_empty in fq_list if non_empty)
    fq_list = sorted(list(itertools.chain(*fq_list)))
    if len(fq_list) == 0:
        print("Allowed R1 patterns:")
        for pattern in read_pattern_list:
            print(pattern)
        raise Exception(
            "\n"
            f"Invalid Read{read} path! \n"
            f"library_id: {library_id}\n"
            f"library_path: {library_path}\n"
        )
    return fq_list


def get_fq(library_id, library_path):
    fq1_list = get_read(library_id, library_path, read="1")
    fq2_list = get_read(library_id, library_path, read="2")
    if len(fq1_list) != len(fq2_list):
        raise Exception("Read1 and Read2 fastq number do not match!")
    fq1 = ",".join(fq1_list)
    fq2 = ",".join(fq2_list)
    return fq1, fq2
