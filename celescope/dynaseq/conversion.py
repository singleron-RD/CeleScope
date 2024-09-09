import pysam
import os
import subprocess
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool

from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.rna.mkref import Mkref_rna
from celescope.tools import reference
from celescope.__init__ import HELP_DICT
from celescope.tools.plotly_plot import Conversion_plot


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
        self.snp_min_cells = args.snp_min_cells
        self.snp_min_depth = args.snp_min_depth
        self.cellsplit = args.cellsplit

        # set
        gtf_file = Mkref_rna.get_config(args.genomeDir)["files"]["gtf"]
        gp = reference.GtfParser(gtf_file)
        gp.get_id_name()
        self.strand_dict = gp.get_strand()
        self.cell_dict, self.cell_num = utils.barcode_list_stamp(
            self.bcfile, cut=self.cellsplit
        )
        self.bam_list = []
        self.conv_df = pd.DataFrame()
        self.snp_df = pd.DataFrame()

        # output files
        ## tmp outputdir
        self.tmp_dir = f"{args.outdir}/tmp"
        utils.check_mkdir(self.tmp_dir)

        ## final outputs
        self.outfile_bam = os.path.join(args.outdir, args.sample + ".PosTag.bam")
        self.outtmp_bam = os.path.join(args.outdir, args.sample + ".PosTag.bam.tmp.bam")
        self.outfile_csv = os.path.join(args.outdir, args.sample + ".PosTag.csv")
        self.outsnp_csv = os.path.join(args.outdir, args.sample + ".snp.csv")

    @utils.add_log
    def run(self):
        # Adding tags and parse snps
        dfs = self.run_conversion()
        # Obtaining conversion positions
        self.snp_candidate(dfs)
        # merge bam files
        self.output_bam()
        # stat and plot
        self.add_conversion_metrics()
        self.conversion_plot()
        # delete tmp dir
        self.clean_tmp()

    def run_cmd(self, cmd):
        subprocess.call(" ".join(cmd), shell=True)

    @utils.add_log
    def run_conversion(self):
        cell_arr, bam_arr = [], []
        fetch_arr = [self.inbam] * len(self.cell_dict)
        strand_arr = [self.strand_dict] * len(self.cell_dict)
        qual_arr = [self.qual] * len(self.cell_dict)
        for x in self.cell_dict:
            tmpbamfile = f"{self.tmp_dir}/tmp_{x}.bam"
            bam_arr.append(tmpbamfile)
            cell_arr.append(self.cell_dict[x])
        self.bam_list = bam_arr

        mincpu = min(self.cell_num, self.thread)
        with Pool(mincpu) as pool:
            results = pool.starmap(
                Conversion.addTags,
                zip(fetch_arr, bam_arr, cell_arr, strand_arr, qual_arr),
            )
        return results

    @utils.add_log
    def snp_candidate(self, df_arr):
        Outputdf = pd.DataFrame()
        for i in df_arr:
            Outputdf = pd.concat([Outputdf, i])
        Outputdf = Outputdf.reset_index()
        # all conv sites
        self.conv_df = Outputdf.groupby("index").agg({"convs": "sum", "cells": "sum"})
        self.conv_df = self.conv_df.reset_index()
        self.conv_df[["chrom", "pos"]] = self.conv_df["index"].str.split(
            "+", expand=True
        )
        self.conv_df = self.conv_df.set_index(["chrom", "pos"])
        self.conv_df.drop("index", axis=1, inplace=True)
        # snp sites
        if self.snp_min_cells < 1:
            self.snp_min_cells = int(self.snp_min_cells * self.cell_num)
        self.snp_df = self.conv_df[self.conv_df["cells"] >= self.snp_min_cells]
        self.snp_df = self.snp_df[self.snp_df["convs"] >= self.snp_min_depth]
        # output
        self.conv_df.to_csv(self.outfile_csv)
        self.snp_df.to_csv(self.outsnp_csv)

    @utils.add_log
    def output_bam(self):
        if len(self.bam_list) > 1:
            bam_list = " ".join(self.bam_list)
            cmd = [
                f"samtools merge -@ {self.thread} -o {self.outtmp_bam}",
                f"{bam_list}",
            ]
            self.run_cmd(cmd)
        else:
            self.outtmp_bam = self.bam_list[0]
        cmd = [
            f"samtools sort -@ {self.thread} -o {self.outfile_bam}",
            f"{self.outtmp_bam}",
        ]
        self.run_cmd(cmd)
        cmd = ["rm", f"{self.outtmp_bam}"]
        self.run_cmd(cmd)

    @utils.add_log
    def clean_tmp(self):
        cmd = f"rm -rf {self.tmp_dir}"
        self.debug_subprocess_call(cmd)

    @staticmethod
    def createTag(d):
        return "".join(["".join(key) + str(d[key]) + ";" for key in d.keys()])[:-1]

    @staticmethod
    def convInRead(read, qual=20):
        specific_conversions = {}
        total_content = {"a": 0, "c": 0, "g": 0, "t": 0}
        specific_conversions[("c", "A")] = 0
        specific_conversions[("g", "A")] = 0
        specific_conversions[("t", "A")] = 0
        specific_conversions[("a", "C")] = 0
        specific_conversions[("g", "C")] = 0
        specific_conversions[("t", "C")] = 0
        specific_conversions[("a", "G")] = 0
        specific_conversions[("c", "G")] = 0
        specific_conversions[("t", "G")] = 0
        specific_conversions[("a", "T")] = 0
        specific_conversions[("c", "T")] = 0
        specific_conversions[("g", "T")] = 0
        specific_conversions[("a", "N")] = 0
        specific_conversions[("c", "N")] = 0
        specific_conversions[("g", "N")] = 0
        specific_conversions[("t", "N")] = 0

        tC_loc = []
        aG_loc = []

        try:
            refseq = read.get_reference_sequence().lower()
        except UnicodeDecodeError:
            return 0
        except AssertionError:
            return 0

        for base in total_content.keys():
            total_content[base] += refseq.count(base)
        for pair in read.get_aligned_pairs(with_seq=True):
            try:
                if pair[0] is not None and pair[1] is not None and pair[2] is not None:
                    if (
                        str(pair[2]).islower()
                        and not read.query_qualities[pair[0]] < qual
                    ):
                        specific_conversions[(pair[2], read.seq[pair[0]])] += 1
                        if (pair[2], read.seq[pair[0]]) == ("t", "C"):
                            tC_loc.append(pair[1])
                        if (pair[2], read.seq[pair[0]]) == ("a", "G"):
                            aG_loc.append(pair[1])
            except (UnicodeDecodeError, KeyError):
                continue
        SC_tag = Conversion.createTag(specific_conversions)
        TC_tag = Conversion.createTag(total_content)

        if len(tC_loc) == 0:
            tC_loc.append(0)
        if len(aG_loc) == 0:
            aG_loc.append(0)
        return SC_tag, TC_tag, tC_loc, aG_loc

    @staticmethod
    def addTags(bamfilename, tmpoutbam, cell_list, strandedness, qual=20):
        tmp_cell = defaultdict(set)
        site_cell = defaultdict(int)
        site_depth = defaultdict(int)
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bamfilename, "rb")
        header = bamfile.header
        mod_bamfile = pysam.AlignmentFile(
            tmpoutbam, mode="wb", header=header, check_sq=False
        )
        pysam.set_verbosity(save)

        class GeneError(Exception):
            pass

        for read in bamfile.fetch(until_eof=True):
            try:
                ## check read info
                if (not read.has_tag("GX")) or read.get_tag("GX") == "-":
                    continue
                if read.get_tag("CB") not in cell_list:
                    continue
                if read.get_tag("GX") not in strandedness:
                    raise GeneError

                tags = Conversion.convInRead(read, qual)
                if tags == 0:
                    continue
                read.set_tag("SC", tags[0], "Z")
                read.set_tag("TC", tags[1], "Z")
                read.set_tag("TL", tags[2])
                read.set_tag("AL", tags[3])
                read.set_tag("ST", strandedness[read.get_tag("GX")])
                mod_bamfile.write(read)

                if strandedness[read.get_tag("GX")] == "+":
                    locs = tags[2]
                else:
                    locs = tags[3]
                if locs[0] != 0:
                    for loc in locs:
                        site = f"{read.reference_name}+{loc}"
                        tmp_cell[site].add(read.get_tag("CB"))
                        site_depth[site] += 1

            except (ValueError, KeyError):
                continue
            except GeneError:
                print(
                    "{} is not in gtf file, please check your files.".format(
                        read.get_tag("GX")
                    )
                )
                continue
            # except (Exception):
            #    print("convert error")
            #    sys.exit(1)
        bamfile.close()
        mod_bamfile.close()

        for x in tmp_cell:
            site_cell[x] = len(tmp_cell[x])
        df1 = pd.DataFrame.from_dict(site_depth, orient="index")
        df2 = pd.DataFrame.from_dict(site_cell, orient="index")
        df = pd.concat([df1, df2], axis=1)
        df.columns = ["convs", "cells"]
        return df

    @utils.add_log
    def add_conversion_metrics(self):
        self.add_metric(
            name="Conversion sites",
            value=self.conv_df.shape[0],
            help_info="the number of T_to_C conversion sites",
        )
        self.add_metric(
            name="Possible SNP sites",
            value=self.snp_df.shape[0],
            help_info="the number of possible SNP sites according to the conditions",
        )

    @utils.add_log
    def conversion_plot(self):
        self.conv_df["cell_pct"] = self.conv_df["cells"] / self.cell_num
        subplot = Conversion_plot(df_bar=self.conv_df).get_plotly_div()
        self.add_data(conversion_box=subplot)


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
        "--snp_min_cells",
        default=10,
        type=float,
        help="Minimum number of cells to call a variant(>=1 for cell number or <1 for cell fraction), default 10.",
    )
    parser.add_argument(
        "--snp_min_depth", default=20, type=int, help="Minimum depth to call a variant"
    )
    parser.add_argument(
        "--cellsplit", default=300, type=int, help="split N cells into a list"
    )
    parser.add_argument(
        "--conversionMem",
        default=30,
        type=int,
        help="Default `30`. Set conversion memory.",
    )

    if sub_program:
        parser.add_argument(
            "--bam",
            help='featureCount bam(sortedByCoord), must have "MD" tag, set in star step',
            required=True,
        )
        parser.add_argument("--cell", help="barcode cell list", required=True)
        parser = s_common(parser)
    return parser
