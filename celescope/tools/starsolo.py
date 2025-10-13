import os
from pathlib import Path
import sys
import subprocess
from collections import Counter

import pandas as pd
import pysam
import celescope.tools.parse_chemistry as parse_chemistry
from celescope.chemistry_dict import chemistry_dict
from typing import Union
from celescope.tools.__init__ import (
    FILTERED_MATRIX_DIR_SUFFIX,
    COUNTS_FILE_NAME,
)
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.tools.barcode import Barcode
from celescope.tools import utils
from celescope.tools.make_ref import MakeRef
from celescope.tools.matrix import CountMatrix
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools.cells import Cells_metrics

SAM_attributes = "NH HI nM AS CR UR CB UB GX GN "
MIN_CELL = 500
MAX_CELL = 60000
MIN_BEAD = 100000
MAX_BEAD = 300000


def create_pattern_args(pattern_dict: dict) -> dict:
    """Create starsolo args related to pattern and return as a dictionary."""
    if len(pattern_dict["U"]) != 1:
        raise ValueError(
            f"Error: Wrong pattern:{pattern_dict}. \n Solution: fix pattern so that UMI only have 1 position.\n"
        )
    ul = pattern_dict["U"][0].start
    ur = pattern_dict["U"][0].stop
    umi_len = ur - ul

    if len(pattern_dict["C"]) == 1:
        solo_type = "CB_UMI_Simple"
        start, stop = pattern_dict["C"][0].start, pattern_dict["C"][0].stop
        cb_start = start + 1
        cb_len = stop - start
        umi_start = ul + 1
        result = {
            "soloType": solo_type,
            "soloCBstart": cb_start,
            "soloCBlen": cb_len,
            "soloUMIstart": umi_start,
            "soloUMIlen": umi_len,
        }
    else:
        solo_type = "CB_UMI_Complex"
        cb_pos = " ".join([f"0_{x.start}_0_{x.stop-1}" for x in pattern_dict["C"]])
        umi_pos = f"0_{ul}_0_{ur-1}"
        result = {
            "soloType": solo_type,
            "soloCBposition": cb_pos,
            "soloUMIposition": umi_pos,
            "soloUMIlen": umi_len,
        }

    return result


def create_v3_pattern_args() -> dict:
    """Can also be used to create flv_rna-V2
    https://github.com/alexdobin/STAR/issues/1607#issuecomment-1210096785
    2: adapter start
    3: adapter end
    """
    linker1 = "ACGATG"
    linker2 = "CATAGT"
    bc = "N" * 9
    linker_len = 6
    bc2_start = 9 + linker_len
    result = {
        "soloType": "CB_UMI_Complex",
        "soloCBposition": " ".join(
            ["2_0_2_8", f"2_{bc2_start}_2_{bc2_start+8}", "3_1_3_9"]
        ),
        "soloUMIposition": "3_10_3_21",
        "soloAdapterSequence": f"{bc}{linker1}{bc}{linker2}",
        "soloAdapterMismatchesNmax": 1,
    }
    return result


def create_whitelist_args(whitelist_str) -> dict:
    """Create whitelist arguments."""
    if not whitelist_str:
        whitelist = "None"
    else:
        whitelist = whitelist_str.strip()
        # nextflow copy remote file to current folder, so only keep file name.
        if whitelist.startswith("http"):
            whitelist = whitelist.split("/")[-1]
        if whitelist.endswith(".gz"):
            whitelist = f"<(gzip -cdf {whitelist})"
    result = {"soloCBwhitelist": whitelist}
    return result


def create_solo_str(
    pattern_args: dict,
    whitelist_args: dict,
    outFileNamePrefix: str,
    fq1: str,
    fq2: str,
    genomeDir: str,
    soloCellFilter: str,
    runThreadN: Union[str, int],
    clip3pAdapterSeq: str,
    outFilterMatchNmin: Union[str, int],
    soloFeatures: str,
    outSAMtype: str,
    outSAMattributes: str,
    soloCBmatchWLtype: str,
    limitBAMsortRAM: Union[str, int],
    extra_starsolo_args: str,
) -> str:
    """
    Create all starsolo args as a dictionary.
    extra_starsolo_args overrides any default args.
    """
    read_command = "zcat" if fq1.strip().endswith(".gz") else "cat"
    args_dict = {
        "outFileNamePrefix": outFileNamePrefix,
        "readFilesIn": f"{fq2} {fq1}",
        "readFilesCommand": read_command,
        "genomeDir": genomeDir,
        "soloCellFilter": soloCellFilter,
        "runThreadN": runThreadN,
        "clip3pAdapterSeq": clip3pAdapterSeq,
        "outFilterMatchNmin": outFilterMatchNmin,
        "soloFeatures": soloFeatures,
        "soloCBmatchWLtype": soloCBmatchWLtype,
        "limitBAMsortRAM": limitBAMsortRAM,
        "outSAMtype": outSAMtype,
        "soloCellReadStats": "Standard",
        "soloBarcodeReadLength": 0,
    }
    for arg in [pattern_args, whitelist_args]:
        args_dict.update(arg)
    if outSAMtype != "None":
        args_dict["outSAMattributes"] = outSAMattributes
    if args_dict["soloType"] == "CB_UMI_Simple":
        args_dict["soloCBmatchWLtype"] = "1MM"
    args_dict.update(cmd_to_dict(extra_starsolo_args))
    cmd = "STAR \\\n"
    for key, value in args_dict.items():
        cmd += f"--{key} {value} \\\n"
    return cmd


def cmd_to_dict(cmd: str) -> dict:
    res = {}
    for part in cmd.split("--"):
        if not part.strip():
            continue
        attr = part.strip().split()
        key = attr[0]
        value = " ".join(attr[1:]) if len(attr) > 1 else ""
        res[key] = value
    return res


def create_soloFeatures(args_soloFeatures, report_soloFeature):
    """
    BAM uses first feature in soloFeatures, so move report_soloFeature to the first position.
    GeneFull_Ex50pAS is also need for intron region metrics.
    """
    features = args_soloFeatures.split(" ")
    if report_soloFeature in features:
        features.remove(report_soloFeature)
    features = [report_soloFeature] + features
    if "GeneFull_Ex50pAS" not in features:
        features.append("GeneFull_Ex50pAS")
    return " ".join(features)


class Starsolo(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        if len(self.fq1_list) != len(self.fq2_list):
            sys.exit("fq1 and fq2 must have same number of files")

        self.chemistry = parse_chemistry.get_chemistry(
            self.assay, args.chemistry, self.fq1_list
        )
        self.pattern_dict, self.bc = parse_chemistry.get_pattern_dict_and_bc(
            self.chemistry, args.pattern, args.whitelist
        )

        self.pattern_args = (
            create_pattern_args(self.pattern_dict)
            if self.chemistry not in ("GEXSCOPE-V3", "flv_rna-V2")
            else create_v3_pattern_args()
        )
        whitelist_str = " ".join(self.bc)
        self.whitelist_args = create_whitelist_args(whitelist_str)
        self.outSAMattributes = SAM_attributes + self.args.SAM_attributes
        self.extra_starsolo_args = args.STAR_param

        # output files
        self.solo_out_dir = f"{self.outdir}/{self.sample}_Solo.out/"
        solo_dir = f"{self.outdir}/{self.sample}_Solo.out/{args.report_soloFeature}"
        self.raw_matrix = Path(f"{solo_dir}/raw")
        self.filtered_matrix = Path(f"{solo_dir}/filtered")
        self.summary_file = f"{solo_dir}/Summary.csv"
        bam = Path(f"{self.outdir}/{self.sample}_Aligned.sortedByCoord.out.bam")

        # outs
        self.outs = [self.raw_matrix, self.filtered_matrix, bam]

    def run_starsolo(self):
        soloFeatures = create_soloFeatures(
            self.args.soloFeatures, self.args.report_soloFeature
        )
        cmd = create_solo_str(
            pattern_args=self.pattern_args,
            whitelist_args=self.whitelist_args,
            outFileNamePrefix=self.out_prefix + "_",
            fq1=self.args.fq1,
            fq2=self.args.fq2,
            genomeDir=self.args.genomeDir,
            soloCellFilter=self.args.soloCellFilter,
            runThreadN=self.args.thread,
            clip3pAdapterSeq=self.args.adapter_3p,
            outFilterMatchNmin=self.args.outFilterMatchNmin,
            soloFeatures=soloFeatures,
            outSAMtype=self.args.outSAMtype,
            outSAMattributes=self.outSAMattributes,
            soloCBmatchWLtype=self.args.soloCBmatchWLtype,
            limitBAMsortRAM=self.args.limitBAMsortRAM,
            extra_starsolo_args=self.extra_starsolo_args,
        )
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)
        cmd = f"chmod -R 755 {self.solo_out_dir}"
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def gzip_matrix(self):
        if not os.path.exists(self.filtered_matrix):
            sys.exit(
                f"{self.filtered_matrix} not found. This error indicates that StarSolo did not generate the gene expression matrix correctly. The specific cause of this error can usually be found in the log above. A common cause of the error is that the fastq file was truncated due to incomplete download."
            )
        cmd = f"gzip {self.raw_matrix}/*; gzip {self.filtered_matrix}/*"
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_Q30_cb_UMI(self):
        fq1_list = self.args.fq1.split(",")
        pattern_dict = self.pattern_dict
        cb_10k, umi_10k, cb_10k_1000k, umi_10k_1000k = (
            Counter(),
            Counter(),
            Counter(),
            Counter(),
        )
        for fq1_file in fq1_list:
            n = 0
            with pysam.FastxFile(fq1_file, persist=False) as fq1:
                for entry in fq1:
                    n += 1
                    if n > 10**6:
                        break
                    qual: str = entry.quality
                    cb_qual = "".join([qual[slice] for slice in pattern_dict["C"]])
                    umi_qual = "".join([qual[slice] for slice in pattern_dict["U"]])
                    if n <= 10**4:
                        cb_10k.update(cb_qual)
                        umi_10k.update(umi_qual)
                    else:
                        cb_10k_1000k.update(cb_qual)
                        umi_10k_1000k.update(umi_qual)

        cb_qual_counter = cb_10k
        umi_qual_counter = umi_10k
        if cb_10k_1000k:
            cb_qual_counter = cb_10k_1000k
            umi_qual_counter = umi_10k_1000k
        q30_cb = sum(
            [cb_qual_counter[k] for k in cb_qual_counter if Barcode.chr_to_int(k) >= 30]
        ) / float(sum(cb_qual_counter.values()))
        q30_umi = sum(
            [
                umi_qual_counter[k]
                for k in umi_qual_counter
                if Barcode.chr_to_int(k) >= 30
            ]
        ) / float(sum(umi_qual_counter.values()))
        return q30_cb, q30_umi

    def run(self):
        self.run_starsolo()
        self.gzip_matrix()
        q30_cb, q30_umi = self.get_Q30_cb_UMI()
        return q30_cb, q30_umi, self.chemistry


def starsolo(args):
    with Starsolo(args) as runner:
        q30_cb, q30_umi, chemistry = runner.run()

    with Mapping(args) as runner:
        valid_reads, corrected = runner.run()

    with Cells(args) as runner:
        n_reads, q30_RNA = runner.run(chemistry, valid_reads)

    with Demultiplexing(args) as runner:
        runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)


class Mapping(Step):
    # only add metrics
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        solo_dir = f"{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS"
        self.cellReadsStats = f"{solo_dir}/CellReads.stats"
        self.filtered_matrix = f"{self.outs_dir}/{FILTERED_MATRIX_DIR_SUFFIX}"
        self.counts_file = f"{solo_dir}/{COUNTS_FILE_NAME }"
        self.genome = MakeRef.get_config(args.genomeDir)["meta"]["genome_name"]

        self.outs = [self.counts_file]

    @utils.add_log
    def parse_cellReadsStats(self) -> dict[str, int]:
        df = pd.read_csv(self.cellReadsStats, sep="\t", header=0, index_col=0)
        df = df.iloc[1:,]  # skip first line cb not pass whitelist
        df_count = df.loc[:, ["nUMIunique", "countedU"]]  # keep dataframe format
        df_count.rename(columns={"nUMIunique": "UMI"}, inplace=True)
        df_count.sort_values(by="UMI", ascending=False, inplace=True)
        cbs = CountMatrix.read_barcodes(self.filtered_matrix)
        df_count["mark"] = "UB"
        for cb in cbs:
            df_count.loc[cb, "mark"] = "CB"
        df_count.to_csv(self.counts_file, sep="\t", index=True)

        df = df.loc[
            :,
            [
                "cbMatch",
                "cbPerfect",
                "genomeU",
                "genomeM",
                "exonic",
                "intronic",
                "exonicAS",
                "intronicAS",
                "countedU",
            ],
        ]
        metrics = df.sum().to_dict()
        return metrics

    @utils.add_log
    def add_metrics_to_report(self, metrics: dict[str, int]):
        valid = metrics["cbMatch"]
        perfect = metrics["cbPerfect"]
        corrected = valid - perfect
        genomeU = metrics["genomeU"]
        genomeM = metrics["genomeM"]
        mapped = genomeU + genomeM
        exonic = metrics["exonic"]
        intronic = metrics["intronic"]
        antisense = metrics["exonicAS"] + metrics["intronicAS"]
        intergenic = mapped - exonic - intronic - antisense
        countedU = metrics["countedU"]

        self.add_metric(
            name="genome",
            value=self.genome,
        )

        self.add_metric(
            name="Reads mapped to unique loci",
            value=genomeU / valid,
            value_type="fraction",
            help_info="Reads that mapped uniquely to the genome.",
        )

        self.add_metric(
            name="Reads mapped to multiple loci",
            value=genomeM / valid,
            value_type="fraction",
            help_info="Reads that mapped to multiple loci in the genome",
        )
        unique_transcriptome = countedU / valid
        self.add_metric(
            name="Reads mapped uniquely to Transcriptome",
            value=unique_transcriptome,
            value_type="fraction",
            help_info="Reads that mapped to a unique gene in the transcriptome. These reads are used for UMI counting.",
        )
        self.add_metric(
            name="Mapped Reads assigned to exonic regions",
            value=exonic / mapped,
            value_type="fraction",
            help_info="Reads that assigned to exonic regions of genes",
        )
        self.add_metric(
            name="Mapped Reads assigned to intronic regions",
            value=intronic / mapped,
            value_type="fraction",
            help_info="Reads that assigned to intronic regions of genes",
        )
        self.add_metric(
            name="Mapped Reads assigned to intergenic regions",
            value=intergenic / mapped,
            value_type="fraction",
            help_info="Reads that can not be assigned to a gene will be considered as intergenic reads.",
        )
        self.add_metric(
            name="Mapped Reads assigned Antisense to gene",
            value=antisense / mapped,
            value_type="fraction",
            help_info="Reads that assigned to the opposite strand of genes",
        )
        return valid, corrected

    @utils.add_log
    def run(self):
        metrics = self.parse_cellReadsStats()
        return self.add_metrics_to_report(metrics)


class Cells(Cells_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        solo_dir = f"{self.outdir}/{self.sample}_Solo.out/{args.report_soloFeature}"
        self.summary_file = f"{solo_dir}/Summary.csv"
        self.counts_file = f"{self.outs_dir}/{COUNTS_FILE_NAME}"

    @utils.add_log
    def parse_summary_add_metrics(self, valid_reads):
        df = pd.read_csv(self.summary_file, index_col=0, header=None)
        s = df.iloc[:, 0]
        n_cells = int(s["Estimated Number of Cells"])
        fraction_reads_in_cells = float(s["Fraction of Unique Reads in Cells"])
        mean_used_reads_per_cell = int(s["Mean Reads per Cell"])
        median_umi_per_cell = int(s["Median UMI per Cell"])
        median_genes_per_cell = int(
            s[f"Median {self.args.report_soloFeature} per Cell"]
        )
        total_genes = int(s[f"Total {self.args.report_soloFeature} Detected"])
        saturation = float(s["Sequencing Saturation"])
        n_reads = int(s["Number of Reads"])
        q30_RNA = float(s["Q30 Bases in RNA read"])

        self.add_cells_metrics(
            n_cells,
            fraction_reads_in_cells,
            mean_used_reads_per_cell,
            median_umi_per_cell,
            total_genes,
            median_genes_per_cell,
            saturation,
            valid_reads,
        )
        return n_reads, q30_RNA, fraction_reads_in_cells

    def run(self, chemistry, valid_reads):
        n_reads, q30_RNA, fraction_reads_in_cells = self.parse_summary_add_metrics(
            valid_reads
        )
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.counts_file))
        return n_reads, q30_RNA


class Demultiplexing(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    def run(self, valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA):
        self.add_metric(
            name="Raw Reads", value=n_reads, help_info="total reads from FASTQ files"
        )
        self.add_metric(
            name="Valid Reads",
            value=valid_reads / n_reads,
            value_type="fraction",
            help_info="fraction of reads with valid barcode and UMI",
        )

        self.add_metric(
            name="Corrected Barcodes",
            value=corrected / valid_reads,
            value_type="fraction",
            help_info="fraction of valid reads with corrected barcodes. Barcodes are corrected to the whitelist sequence that is 1 Hamming-distance away.",
        )

        self.add_metric(
            name="Q30 of Barcodes",
            value=q30_cb,
            value_type="fraction",
            help_info="Fraction of barcode bases with quality score >= 30",
        )

        self.add_metric(
            name="Q30 of UMI",
            value=q30_umi,
            value_type="fraction",
            help_info="Fraction of UMI bases with quality score >= 30",
        )

        self.add_metric(
            name="Q30 of RNA Reads",
            value=q30_RNA,
            value_type="fraction",
            help_info="Fraction of RNA read bases with quality score >= 30",
        )


def get_opts_starsolo(parser, sub_program=True):
    parser.add_argument(
        "--chemistry",
        help=HELP_DICT["chemistry"],
        choices=list(chemistry_dict.keys()),
        default="auto",
    )
    parser.add_argument(
        "--pattern",
        help="""The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  
        - `C`: cell barcode  
        - `L`: linker(common sequences)  
        - `U`: UMI    
        - `T`: poly T""",
    )
    parser.add_argument(
        "--whitelist",
        help="Cell barcode whitelist file path, one cell barcode per line. Multiple whitelist are seperated by whitespace.",
    )
    parser.add_argument(
        "--adapter_3p",
        help="Adapter sequence to clip from 3 prime. Multiple sequences are seperated by space",
        default="AAAAAAAAAAAA",
    )
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
        "--soloCellFilter",
        help="Same as the argument in STARsolo",
        default="EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000",
    )
    parser.add_argument(
        "--limitBAMsortRAM",
        help="Same as the argument in STARsolo",
        default=32000000000,
        type=int,
    )
    parser.add_argument("--STAR_param", help=HELP_DICT["additional_param"], default="")
    parser.add_argument(
        "--outSAMtype",
        help="Same as the argument in STARsolo. Set to 'None' to skip BAM file output.",
        default="BAM SortedByCoordinate",
        type=str,
    )
    parser.add_argument(
        "--SAM_attributes",
        help=f"Additional attributes(other than {SAM_attributes}) to be added to SAM file",
        default="",
    )
    parser.add_argument(
        "--soloFeatures",
        help="Same as the argument in STARsolo",
        default="GeneFull_Ex50pAS Gene",
    )
    parser.add_argument(
        "--soloCBmatchWLtype",
        help="Same as the argument in STARsolo. Please note `EditDist_2` only works with `--soloType CB_UMI_Complex`. ",
        default="EditDist_2",
        type=str,
    )
    parser.add_argument(
        "--report_soloFeature",
        help="Specify which soloFeatures to use in the HTML report and the outs directory.",
        default="GeneFull_Ex50pAS",
    )

    if sub_program:
        parser.add_argument(
            "--fq1",
            help="R1 fastq file. Multiple files are separated by comma.",
            required=True,
        )
        parser.add_argument(
            "--fq2",
            help="R2 fastq file. Multiple files are separated by comma.",
        )
        parser = s_common(parser)

    return parser
