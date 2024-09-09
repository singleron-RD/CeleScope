import subprocess
import os

from celescope.tools import utils


class Mkref:
    """
    ## Features

    - Build Index for IMGT_ref.

    - Make sure current directory contains all V,D,J of TRA/TRB or IGH/IGK/IGL reference downloaded from IMGT website.

    ## Usage
    ```
    celescope vdj mkref species(human or mouse) imgt_type(TR or IG)
    ```

    ## Output

    - VDJ IMGT reference with index files.

    """

    def __init__(self, args):
        self.type = args.imgt_type
        self.species = args.species
        self.outdir = f"{os.getcwd()}/{self.species}_{self.type}"

    def __call__(self):
        utils.check_mkdir(self.outdir)
        self.combine_seq()
        self.build_index()

    @utils.add_log
    def combine_seq(self):
        """
        Combine all V, all D and all J sequences, respectively, into separate files:
        """
        if self.type == "TR":
            cmd_V = f"cat TRAV.fasta TRBV.fasta > {self.outdir}/TRV.fasta"
            cmd_D = f"cat TRBD.fasta > {self.outdir}/TRD.fasta"
            cmd_J = f"cat TRAJ.fasta TRBJ.fasta > {self.outdir}/TRJ.fasta"
        elif self.type == "TRGD":
            cmd_V = f"cat TRDV.fasta TRGV.fasta > {self.outdir}/TRV.fasta"
            cmd_D = f"cat TRDD.fasta > {self.outdir}/TRD.fasta"
            cmd_J = f"cat TRDJ.fasta TRGJ.fasta > {self.outdir}/TRJ.fasta"
        else:
            cmd_V = f"cat IGHV.fasta IGKV.fasta IGLV.fasta > {self.outdir}/IGV.fasta"
            cmd_D = f"cat IGHD.fasta > {self.outdir}/IGD.fasta"
            cmd_J = f"cat IGHJ.fasta IGKJ.fasta IGLJ.fasta > {self.outdir}/IGJ.fasta"

        for cmd in [cmd_V, cmd_D, cmd_J]:
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def build_index(self):
        """
        build imgt index for igblast
        """
        os.chdir(self.outdir)
        for file in os.listdir():
            out_file_name, _ = os.path.splitext(file)
            subprocess.check_call(
                f"edit_imgt_file.pl {file} > {out_file_name}.fa", shell=True
            )
            subprocess.check_call(
                f"makeblastdb -parse_seqids -dbtype nucl -in {out_file_name}.fa",
                shell=True,
            )


def mkref(args):
    runner = Mkref(args)
    runner()


def get_opts_mkref(parser, sub_program=True):
    if sub_program:
        parser.add_argument("species", help="human or mouse.")
        parser.add_argument(
            "imgt_type",
            help="TR: build index for TR IMGT reference, IG: build index for IG IMGT reference",
        )

    return parser
