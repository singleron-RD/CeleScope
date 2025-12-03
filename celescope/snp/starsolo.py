import pandas as pd
from celescope.tools.starsolo import (
    Starsolo as tools_Starsolo,
    get_opts_starsolo as tools_opts,
    Mapping as tools_Mapping,
    Demultiplexing,
)
from celescope.tools.step import Step


class Starsolo(tools_Starsolo):
    def __init__(self, args):
        super().__init__(args)
        self.outs = []


class Mapping(tools_Mapping):
    def __init__(self, args):
        super().__init__(args)
        self.outs = []
        solo_dir = f"{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS"
        self.filtered_matrix = f"{solo_dir}/filtered"


class Summary(Step):
    def __init__(self, args):
        super().__init__(args)
        solo_dir = f"{self.outdir}/{self.sample}_Solo.out/{args.report_soloFeature}"
        self.summary_file = f"{solo_dir}/Summary.csv"

    def run(self):
        df = pd.read_csv(self.summary_file, index_col=0, header=None)
        s = df.iloc[:, 0]
        q30_RNA = float(s["Q30 Bases in RNA read"])
        n_reads = int(s["Number of Reads"])
        return n_reads, q30_RNA


def get_opts_starsolo(parser, sub_program=True):
    tools_opts(parser, sub_program)
    return parser


def starsolo(args):
    with Starsolo(args) as runner:
        q30_cb, q30_umi, chemistry = runner.run()

    with Mapping(args) as runner:
        valid_reads, corrected = runner.run()

    with Summary(args) as runner:
        n_reads, q30_RNA = runner.run()

    with Demultiplexing(args) as runner:
        runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)
