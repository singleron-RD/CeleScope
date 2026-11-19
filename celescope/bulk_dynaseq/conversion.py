import os
import glob
import pandas as pd
from multiprocessing import Pool

from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.rna.mkref import Mkref_rna
from celescope.tools import reference
from celescope.__init__ import HELP_DICT
from celescope.dynaseq.conversion_worker import conversion_process


class Conversion(Step):
    """
    ## Features
    - Get conversion pos in each read using well-split bams from the split step.
        - Get snp info.
    - Per-well outputs.

    ## Output
    - `{well}/{well}.PosTag.bam` Per-well bam file with conversion info.
    - `{well}/{well}.PosTag.csv` Per-well TC conversion sites info in csv format.
    - `{well}/{well}.snp.csv` Per-well candidated snp sites.
    """

    def __init__(self, args):
        Step.__init__(self, args)
        self.sample = args.sample
        self.bam_dir = args.bam_dir  # directory containing well-split bams
        self.bcfile = args.barcodes  # barcode file
        self.outdir = args.outdir
        self.thread = int(args.thread)
        self.qual = int(args.basequalilty)
        self.snp_min_depth = args.snp_min_depth
        self.snp_threshold = args.snp_threshold

        gtf_file = Mkref_rna.get_config(args.genomeDir)["files"]["gtf"]
        gp = reference.GtfParser(gtf_file)
        gp.get_id_name()
        self.strand_dict = gp.get_strand()
        self.cell_list, self.cell_num = utils.read_one_col(self.bcfile)

    @utils.add_log
    def run(self):
        fetch_arr = sorted(glob.glob(f"{self.bam_dir}/*.bam", recursive=False))
        if len(fetch_arr) == 0:
            raise ValueError(
                f"No bam files found in {self.bam_dir}. "
                "Please check the split step output."
            )

        well_names = []
        bam_arr = []
        for x in fetch_arr:
            well_name = os.path.splitext(os.path.basename(x))[0]
            well_names.append(well_name)
            well_outdir = os.path.join(self.outdir, well_name)
            utils.check_mkdir(well_outdir)
            tmpbamfile = os.path.join(well_outdir, well_name + ".PosTag.bam")
            bam_arr.append(tmpbamfile)

        cell_arr = [self.cell_list] * len(fetch_arr)
        strand_arr = [self.strand_dict] * len(fetch_arr)
        qual_arr = [self.qual] * len(fetch_arr)

        mincpu = min(self.cell_num, self.thread)
        with Pool(mincpu) as pool:
            results = pool.map(
                conversion_process,
                zip(fetch_arr, bam_arr, cell_arr, strand_arr, qual_arr),
            )

        # Per-well output
        for well_name, df in zip(well_names, results):
            self.per_well_output(well_name, df)

    def per_well_output(self, well_name, df):
        well_outdir = os.path.join(self.outdir, well_name)
        if df.empty:
            conv_df = pd.DataFrame(
                columns=["chrom", "pos", "convs", "covers", "posratio"]
            )
            snp_df = pd.DataFrame(
                columns=["chrom", "pos", "convs", "covers", "posratio"]
            )
        else:
            Outputdf = df.reset_index()
            conv_df = Outputdf.groupby(["chrom", "pos"]).agg(
                {"convs": "sum", "covers": "sum"}
            )
            conv_df["posratio"] = conv_df["convs"] / conv_df["covers"]
            conv_df = conv_df.reset_index()
            snp_df = conv_df[conv_df["posratio"] >= self.snp_threshold]
            snp_df = snp_df[snp_df["covers"] >= self.snp_min_depth]

        conv_df.to_csv(
            os.path.join(well_outdir, f"{well_name}.PosTag.csv"), index=False
        )
        snp_df.to_csv(os.path.join(well_outdir, f"{well_name}.snp.csv"), index=False)


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
        "--conversionMem",
        default=5,
        type=int,
        help="Set conversion memory.",
    )
    parser.add_argument(
        "--well_sample",
        help="tsv file of well numbers and sample names. The first column is well numbers and the second column is sample names.",
        required=True,
    )

    if sub_program:
        parser.add_argument(
            "--bam_dir",
            help="Directory containing well-split bam files from the split step",
            required=True,
        )
        parser.add_argument(
            "--barcodes",
            help="barcode file.",
            required=True,
        )
        parser = s_common(parser)
    return parser
