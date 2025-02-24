import subprocess
import sys
from celescope.tools.starsolo import (
    Starsolo as tools_Starsolo,
    create_solo_args,
    Mapping,
    Cells,
    Demultiplexing,
)
from celescope.tools.step import s_common
from celescope.__init__ import HELP_DICT
from celescope.chemistry_dict import chemistry_dict

SAM_attributes = "NH HI nM AS CR UR CB UB GX GN "


class Starsolo(tools_Starsolo):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    def run_starsolo(self):
        cmd = create_solo_args(
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
            soloFeatures=self.args.soloFeatures,
            outSAMattributes=self.outSAMattributes,
            soloCBmatchWLtype=self.args.soloCBmatchWLtype,
            extra_starsolo_args=self.extra_starsolo_args,
        )
        if self.chemistry == "bulk_rna-bulk_vdj_match":
            cmd += "--soloStrand Reverse \\\n"
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)
        cmd = f"chmod -R 755 {self.solo_out_dir}"
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)


def starsolo(args):
    with Starsolo(args) as runner:
        q30_cb, q30_umi, chemistry = runner.run()

    with Mapping(args) as runner:
        valid_reads, corrected = runner.run()

    with Cells(args) as runner:
        n_reads, q30_RNA = runner.run(chemistry, valid_reads)

    with Demultiplexing(args) as runner:
        runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)


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
        help="Cell barcode whitelist file path, one cell barcode per line.",
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
    # 96 wells
    parser.add_argument(
        "--soloCellFilter",
        help="The same as the argument in STARsolo",
        default="CellRanger2.2 384 0.99 10",
    )
    parser.add_argument(
        "--starMem", help="Maximum memory that STAR can use.", default=32
    )
    parser.add_argument("--STAR_param", help=HELP_DICT["additional_param"], default="")
    parser.add_argument(
        "--SAM_attributes",
        help=f"Additional attributes(other than {SAM_attributes}) to be added to SAM file",
        default="",
    )
    parser.add_argument(
        "--soloFeatures",
        help="The same as the argument in STARsolo",
        default="GeneFull_Ex50pAS Gene",
    )
    parser.add_argument(
        "--soloCBmatchWLtype",
        help="The same as the argument in STARsolo. Please note `EditDist 2` only works with `--soloType CB UMI Complex`. ",
        default="1MM",
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
            required=True,
        )
        parser = s_common(parser)

    return parser
