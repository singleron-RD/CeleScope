import os
import glob
import subprocess
import pandas as pd
from multiprocessing import Pool

from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.rna.mkref import Mkref_rna
from celescope.tools import reference
from celescope.__init__ import HELP_DICT
from celescope.dynaseq.conversion_worker import conversion_process


def parse_star_log(star_log_file):
    """
    parse Number of input reads from STAR log

    Args:
        star_log_file (str): STAR log dir

    Returns:
        int: Number of input reads, None if not found
    """
    try:
        with open(star_log_file, "r") as f:
            for line in f:
                if "Number of input reads |" in line:
                    num_reads = line.strip().split("|")[1].strip()
                    return int(num_reads)
    except Exception as e:
        print(f"Error when parsing STAR log file: {e}")
        return None


class Conversion(Step):
    """
    ## Features
    - Get conversion pos in each read.
        - Get snp info.

    ## Output
    - `{sample}.PosTag.bam` Bam file with conversion info.
    - `{sample}.PosTag.csv` TC conversion sites info in csv format.
    - `{sample}.snp.csv` Candidated snp sites.
    """

    def __init__(self, args):
        Step.__init__(self, args)
        # input files
        self.sample = args.sample
        self.inbam = args.bam
        self.bcfile = args.cell
        self.outdir = args.outdir
        self.thread = int(args.thread)
        self.qual = int(args.basequalilty)
        self.snp_min_depth = args.snp_min_depth
        self.snp_threshold = args.snp_threshold
        self.readsplit = args.readsplit

        # set
        gtf_file = Mkref_rna.get_config(args.genomeDir)["files"]["gtf"]
        gp = reference.GtfParser(gtf_file)
        gp.get_id_name()
        self.strand_dict = gp.get_strand()
        self.cell_list, self.cell_num = utils.read_one_col(self.bcfile)
        if args.total_num:
            self.total_num = args.total_num
        elif args.star_log:
            self.total_num = parse_star_log(args.star_log)
        else:
            self.total_num = None

        self.bam_list = []
        self.conv_df = pd.DataFrame()
        self.snp_df = pd.DataFrame()

        # output files
        ## tmp outputdir
        self.tmp_dir = f"{args.outdir}/../tmp/"
        utils.check_mkdir(self.tmp_dir)

        ## final outputs
        self.outfile_bam = os.path.join(args.outdir, args.sample + ".PosTag.bam")
        self.outfile_csv = os.path.join(args.outdir, args.sample + ".PosTag.csv")
        self.outsnp_csv = os.path.join(args.outdir, args.sample + ".snp.csv")

    @utils.add_log
    def run(self):
        # Split bams, CalledProcessError
        self.split_bam()
        # Adding tags and parse snps
        dfs = self.run_conversion()
        # Obtaining conversion positions
        self.snp_candidate(dfs)
        # merge bam files
        self.output_bam()
        # stat and plot
        self.add_conversion_metrics()
        # delete tmp dir
        self.clean_tmp()

    @utils.add_log
    def split_bam(self):
        if self.total_num:
            param_total = f"TOTAL_READS_IN_INPUT={self.total_num}"
        else:
            param_total = ""
        tmp_dir = f"{self.tmp_dir}/0"
        utils.check_mkdir(tmp_dir)
        cmd = (
            "picard SplitSamByNumberOfReads \\\n"
            f"I={self.inbam} \\\n"
            f"O={tmp_dir} \\\n"
            f"OUT_PREFIX={self.sample} \\\n"
            f"SPLIT_TO_N_READS={self.readsplit} \\\n"
            f"{param_total} \\\n"
        )
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run_conversion(self):
        tmp_dir = f"{self.tmp_dir}/1"
        utils.check_mkdir(tmp_dir)
        bam_arr = []
        fetch_arr = glob.glob(f"{self.tmp_dir}/0/*.bam", recursive=False)
        cell_arr = [self.cell_list] * len(fetch_arr)
        strand_arr = [self.strand_dict] * len(fetch_arr)
        qual_arr = [self.qual] * len(fetch_arr)
        for x in fetch_arr:
            tmpbamfile = os.path.join(tmp_dir, os.path.basename(x))
            bam_arr.append(tmpbamfile)
        self.bam_list = bam_arr

        mincpu = min(self.cell_num, self.thread)
        with Pool(mincpu) as pool:
            results = pool.map(
                conversion_process,
                zip(fetch_arr, bam_arr, cell_arr, strand_arr, qual_arr),
            )
        return results

    @utils.add_log
    def snp_candidate(self, df_arr):
        Outputdf = pd.concat(df_arr)
        Outputdf = Outputdf.reset_index()
        # all conv sites
        self.conv_df = Outputdf.groupby(["chrom", "pos"]).agg(
            {"convs": "sum", "covers": "sum"}
        )
        self.conv_df["posratio"] = self.conv_df["convs"] / self.conv_df["covers"]
        self.conv_df = self.conv_df.reset_index()
        # snp sites
        self.snp_df = self.conv_df[self.conv_df["posratio"] >= self.snp_threshold]
        self.snp_df = self.snp_df[self.snp_df["covers"] >= self.snp_min_depth]
        # output
        self.conv_df.to_csv(self.outfile_csv, index=False)
        self.snp_df.to_csv(self.outsnp_csv, index=False)

    @utils.add_log
    def output_bam(self):
        if len(self.bam_list) > 1:
            bam_list = " ".join(self.bam_list)
            cmd = (
                f"samtools merge -f -@ {self.thread} -o {self.outfile_bam} "
                f"{bam_list}"
            )
            subprocess.check_call(cmd, shell=True)
        else:
            cmd = f"mv {self.bam_list[0]} {self.outfile_bam}"
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def clean_tmp(self):
        cmd = f"rm -rf {self.tmp_dir}/0/"
        self.debug_subprocess_call(cmd)

    @utils.add_log
    def add_conversion_metrics(self):
        self.add_metric(
            name="Conversion sites",
            value=self.conv_df.shape[0],
            help_info="The total number of T_to_C conversion sites",
        )
        self.add_metric(
            name="Possible SNP sites",
            value=self.snp_df.shape[0],
            help_info=(
                f"The total number of SNP sites under the specified conditions: "
                f"<code>--snp_threshold {self.snp_threshold}</code> and "
                f"<code>--snp_min_depth {self.snp_min_depth}</code>"
            ),
        )


@utils.add_log
def conversion(args):
    with Conversion(args) as runner:
        runner.run()


def get_opts_conversion(parser, sub_program):
    parser.add_argument("--genomeDir", help=HELP_DICT["genomeDir"])
    parser.add_argument(
        "--basequalilty",
        default=20,
        type=int,
        help="min base quality of the read sequence",
        required=False,
    )
    parser.add_argument(
        "--snp_min_depth", default=20, type=int, help="Minimum depth to call a snp"
    )
    parser.add_argument(
        "--snp_threshold",
        default=0.5,
        type=float,
        help="snp threshold filter, greater than snp_threshold will be recognized as snp",
    )
    parser.add_argument(
        "--readsplit",
        default=1000000,
        type=int,
        help="split to have approximately N reads per output file",
    )
    parser.add_argument(
        "--conversionMem",
        default=5,
        type=int,
        help="Set conversion memory.",
    )

    if sub_program:
        parser.add_argument(
            "--bam",
            help='starsolo output bam(sortedByCoord), must have "MD" tag, set in starsolo step',
            required=True,
        )
        parser.add_argument("--cell", help="barcode cell list", required=False)
        parser.add_argument(
            "--total_num", help="total number of input bam", type=int, required=False
        )
        parser.add_argument(
            "--star_log",
            help="star log for extracting the total number (alternative to `total_num`)",
            required=False,
        )
        parser = s_common(parser)
    return parser
