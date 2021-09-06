import subprocess
import os

from celescope.tools import utils
from celescope.trust_vdj.__init__ import TRUST, INDEX
from celescope.tools.step import Step, s_common

class Assemble(Step):
    """
    Features

    - Assemble TCR/BCR seq data.

    Output

    - `03.assemble/{sample}_toassemble.fq` Reads to assemble.
    - `03.assemble/{sample}_toassemble_bc.fa` Barcodes to assemble.
    - `03.assemble/{sample}_cdr3.out` All assembled CDR3 output.
    - `03.assemble/{sample}_barcode_report.tsv` Record chain information in each barcode.
    - `03.assemble/{sample}_annot.fa` Assembled annotated contig sequences.
    - `03.assemble/{sample}_assembled_reads.fa` Assembled raw reads.
    - `03.assemble/{sample}_report.tsv` Record assembled CDR3 types and count.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = f'{args.outdir}/all_outs'
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.sample = args.sample
        self.species = args.species
        self.speed_up = args.speed_up

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @utils.add_log
    def run(self):
        species = self.species
        index_file = f'{INDEX}/{species}/{species}_ref.fa'
        ref = f'{INDEX}/{species}/{species}_IMGT+C.fa'
        string1 = ''
        if self.speed_up:
            string1 = '--repseq '

        cmd = (
            f'run-trust4 -t {self.thread} '
            f'-u {self.fq2} '
            f'--barcode {self.fq1} '
            f'--barcodeRange 0 23 + '
            f'-f {index_file} '
            f'--ref {ref} '
            f'{string1}'
            f'-o {self.sample} --od {self.outdir}' 
        )

        Assemble.run.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_full_len_assembly(self):
        cmd = (
            f'perl {TRUST}/GetFullLengthAssembly.pl {self.outdir}/{self.sample}_annot.fa > '
            f'{self.outdir}/{self.sample}_full_len.fa '
        )
        Assemble.get_full_len_assembly.info(cmd)
        subprocess.check_call(cmd, shell=True)

@utils.add_log
def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='R1 reads matched with scRNA-seq.', required=True)
        parser.add_argument('--fq2', help='R2 reads matched with scRNA-seq.', required=True)

    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)
    parser.add_argument('--speed_up', help='Speed assemble for TCR/BCR seq data.', action='store_true')    








