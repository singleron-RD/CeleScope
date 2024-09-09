import subprocess

import pandas as pd

from celescope.tools import utils
from celescope.tools.make_ref import MakeRef
from celescope.tools.step import s_common
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step
from celescope.tools.__init__ import STAR_BAM_SUFFIX


@utils.add_log
def get_star_cmd(args, input_file, output_prefix):
    """
    output sam format to improve speed
    """
    cmd = (
        f"STAR "
        f"--runThreadN {args.thread} "
        f"--genomeDir {args.genomeDir} "
        f"--outSAMmultNmax 1 "
        f"--outFilterMultimapNmax {args.outFilterMultimapNmax} "
        f"--outSAMtype BAM Unsorted "
        f"--outFilterMatchNmin {args.outFilterMatchNmin} "
        f"{args.STAR_param} "
        f"--readFilesIn {input_file} "
        f"--outFileNamePrefix {output_prefix}_ "
    )
    get_star_cmd.logger.info(cmd)
    return cmd


def get_star_log(logPath):
    df = pd.read_csv(logPath, sep="\t", header=None, names=["name", "value"])
    log_dict = {}
    for t in df.itertuples():
        name = t.name.strip("|").strip()
        log_dict[name] = t.value

    return log_dict


class Star_mixin(Step):
    """
    Mixin class for STAR
    """

    def __init__(self, args, add_prefix=None, display_title=None):
        super().__init__(args, display_title)

        # parse
        self.genome_name = MakeRef.get_config(args.genomeDir)["meta"]["genome_name"]
        self.stat_prefix = "Reads"
        if getattr(args, "consensus_fq", False):
            self.stat_prefix = "UMIs"

        # out
        if add_prefix:
            self.out_prefix += f"_{add_prefix}"
        self.logPath = f"{self.out_prefix}_Log.final.out"
        self.STAR_bam = f"{self.out_prefix}_{STAR_BAM_SUFFIX}"

    @utils.add_log
    def STAR(self):
        input_file = self.args.fq
        output_prefix = self.out_prefix
        cmd = get_star_cmd(self.args, input_file, output_prefix)

        if self.args.out_unmapped:
            cmd += "--outReadsUnmapped Fastx"
        if self.args.fq[-3:] == ".gz":
            cmd += "--readFilesCommand zcat"
        if self.args.STAR_param:
            cmd += " " + self.args.STAR_param
        self.STAR.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run(self):
        self.STAR()
        self.add_star_metrics()

    def add_star_metrics(self):
        """
        step metrics
        """
        log_dict = get_star_log(self.logPath)

        unique_reads = int(log_dict["Uniquely mapped reads number"])
        multiple_loci_reads = int(log_dict["Number of reads mapped to multiple loci"])
        too_many_loci_reads = int(log_dict["Number of reads mapped to too many loci"])
        total_reads = int(log_dict["Number of input reads"])
        multi_reads = multiple_loci_reads + too_many_loci_reads

        self.add_metric(
            name="Genome",
            value=self.genome_name,
        )
        self.add_metric(
            name=f"Uniquely Mapped {self.stat_prefix}",
            value=unique_reads,
            total=total_reads,
            help_info="reads that mapped uniquely to the genome",
        )
        self.add_metric(
            name=f"Multi-Mapped {self.stat_prefix}",
            value=multi_reads,
            total=total_reads,
            help_info="reads that mapped to multiple locations in the genome",
        )


def get_opts_star_mixin(parser, sub_program):
    parser.add_argument(
        "--genomeDir",
        help=HELP_DICT["genomeDir"],
    )
    parser.add_argument(
        "--outFilterMatchNmin",
        help="""Alignment will be output only if the number of matched bases 
is higher than or equal to this value.""",
        default=50,
    )
    parser.add_argument(
        "--out_unmapped", help="Output unmapped reads.", action="store_true"
    )
    parser.add_argument("--STAR_param", help=HELP_DICT["additional_param"], default="")
    parser.add_argument(
        "--outFilterMultimapNmax",
        help=(
            "How many places are allowed to match a read at most. Please note that even if this value is greater than 1, "
            "at most 1 alignment(with the highest score) will be output in the bam file."
        ),
        default=1,
    )
    parser.add_argument(
        "--starMem", help="Default `30`. Maximum memory that STAR can use.", default=30
    )
    if sub_program:
        parser.add_argument("--fq", help="Required. R2 fastq file.", required=True)
        parser.add_argument(
            "--consensus_fq",
            action="store_true",
            help="A indicator that the input fastq has been consensused.",
        )
        parser = s_common(parser)
