import argparse
import sys
import subprocess

from celescope.tools.step import Step
from celescope.tools.barcode import get_opts_barcode
from celescope.tools.cutadapt import get_opts_cutadapt, get_cutadapt_cmd, Cutadapt as Ct
from celescope.rna.star import get_opts_star
from celescope.tools.star_mixin import get_star_cmd, Star_mixin
from celescope.tools import utils


def args_to_list(args):
    s = []
    for arg, v in vars(args).items():
        if v not in (False, None, "") and arg != "func":
            s.append(f"--{arg}")
            if v is not True:
                if "param" in arg:
                    s.append(f'"{v}"')
                else:
                    s.append(f"{v}")
    return s


def get_step_args_str(args, get_opts_func):
    args_str = ""
    parser = argparse.ArgumentParser()
    get_opts_func(parser, sub_program=True)
    # [known_args_namespace, not_known_args_list]
    known_args = parser.parse_known_args(args_to_list(args))[0]
    args_str = " ".join(args_to_list(known_args))
    return args_str


class Prep(Step):
    def __init__(self, args):
        super().__init__(args)

    @utils.add_log
    def run(self):
        bc2 = f"{self.out_prefix}_bc2.FIFO"
        ct2 = f"{self.out_prefix}_ct2.FIFO"
        bc_args_str = get_step_args_str(self.args, get_opts_barcode)
        ct_cmd = get_cutadapt_cmd(self.args, bc2, ct2)
        star_cmd = get_star_cmd(self.args, ct2, self.out_prefix)
        cmd = (
            f"mkfifo {bc2} {ct2} \n"
            f"celescope rna barcode {bc_args_str} --stdout >> {bc2} & \n"
            f"{ct_cmd} & \n"
            f"{star_cmd} \n"
            f"rm {bc2} {ct2} \n"
        )
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)


class Cutadapt(Ct):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    def run(self):
        self.add_cutadapt_metrics()


class Star(Star_mixin):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    def run(self):
        self.sort_bam()
        self.index_bam()
        self.remove_unsort_bam()
        self.add_star_metrics()


@utils.add_log
def prep(args):
    runner = Prep(args)
    runner.run()

    with Cutadapt(args, display_title="Trimming") as runner:
        runner.run()

    with Star(args, display_title="Mapping") as runner:
        runner.run()


def get_opts_prep(parser, sub_program):
    get_opts_barcode(parser, sub_program=sub_program)
    get_opts_cutadapt(parser, sub_program=False)
    get_opts_star(parser, sub_program=False)
