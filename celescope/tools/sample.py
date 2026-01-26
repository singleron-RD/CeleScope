from pathlib import Path
from celescope.tools import utils
from celescope.__init__ import __VERSION__
from celescope.chemistry_dict import chemistry_dict
import subprocess
import tarfile

from celescope.tools.step import Step, s_common

from celescope.tools.parse_chemistry import get_chemistry, invalid_debug


def sub_fq(fq, out_fq, use_read=1000, skip_read=10000):
    cmd = f"zcat {fq} | head -n {(use_read + skip_read) * 4} | tail -n {use_read * 4} | gzip > {out_fq}"
    subprocess.check_call(cmd, shell=True)


class Sample(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        # out
        self.invalid_debug_file = f"{self.out_prefix}_invalid_debug.html"
        self.subsample_fqs: list[Path] = []

    @utils.add_log
    def run(self):
        chemistry = get_chemistry(self.assay, self.args.chemistry, self.fq1_list)

        self.add_metric(
            name="Sample ID",
            value=self.sample,
        )
        self.add_metric(
            name="Assay",
            value=self.assay,
        )
        self.add_metric(
            name="Chemistry",
            value=chemistry,
            help_info='For more information, see <a href="https://github.com/singleron-RD/CeleScope/blob/master/doc/chemistry.md">here</a>',
        )
        self.add_metric(
            name="Software Version",
            value=__VERSION__,
        )
        invalid_debug(
            chemistry,
            self.fq1_list,
            self.invalid_debug_file,
            pattern=self.args.pattern,
            whitelist=self.args.whitelist,
            linker=self.args.linker,
            use_read=self.args.use_read,
            skip_read=self.args.skip_read,
        )
        self.subsample(self.fq1_list)
        self.subsample(self.fq2_list)
        self.pack_files()

    def subsample(self, fq_list):
        for fq in fq_list:
            fq = Path(fq)
            out_fq = Path(self.outdir) / f"{fq.stem}.subsample.gz"
            sub_fq(fq, out_fq, self.args.use_read, self.args.skip_read)
            self.subsample_fqs.append(out_fq)

    def pack_files(self):
        with tarfile.open(f"{self.out_prefix}_debug.tar", "w") as tar:
            tar.add(self.invalid_debug_file, arcname=Path(self.invalid_debug_file).name)
            for file in self.subsample_fqs:
                tar.add(file, arcname=file.name)


@utils.add_log
def sample(args):
    with Sample(args) as runner:
        runner.run()


def get_opts_sample(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument("--fq1", help="read1 fq file")
        parser.add_argument("--fq2", help="read2 fq file")
        parser.add_argument("--linker")
        parser.add_argument(
            "--use_read", help="number of reads used for debug", default=1000, type=int
        )
        parser.add_argument(
            "--skip_read",
            help="Number of reads to skip for debugging. The R1 FASTQ contains N bases in the first hundred reads, so it is better to skip them during debugging",
            default=10000,
            type=int,
        )
    parser.add_argument(
        "--chemistry",
        choices=list(chemistry_dict.keys()),
        help="chemistry version",
        default="auto",
    )
    parser.add_argument("--pattern")
    parser.add_argument(
        "--whitelist",
    )
    return parser
