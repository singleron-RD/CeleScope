from celescope.tools import utils
from celescope.tools.step import Step, s_common
import subprocess

SOFTWARE = '/SGRNJ03/randd/zhouxin/software/TRUST4/'



class Assemble(Step):
    """
    Features

    - Assemble TCR/BCR seq data.

    Output

    - `05.assemble/{sample}_toassemble.fq` Reads to assemble.
    - `05.assemble/{sample}_toassemble_bc.fa` Barcodes to assemble.
    - `05.assemble/{sample}_cdr3.out` All assembled CDR3 output.
    - `05.assemble/{sample}_barcode_report.tsv` Record chain information in each barcode.
    - `05.assemble/{sample}_annot.fa` Assembled annotated contig sequences.
    - `05.assemble/{sample}_assembled_reads.fa` Assembled raw reads.
    - `05.assemble/{sample}_report.tsv` Record assembled CDR3 types and count.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.sample = args.sample
        self.species = args.species
        self.speed_up = args.speed_up


    @utils.add_log
    def run(self):

        species = self.species

        index_file = f'{SOFTWARE}/index/{species}/{species}_ref.fa'
        ref = f'{SOFTWARE}/index/{species}/{species}_IMGT+C.fa'

        string1 = ''
        if self.speed_up:
            string1 = '--repseq '
        cmd = (
            f'source activate zhouxinT; '
            f'{SOFTWARE}/run-trust4 -t {self.thread} '
            f'-u {self.fq2} '
            f'--barcode {self.fq1} '
            f'--barcodeRange 0 23 + '
            f'--UMI {self.fq1} '
            f'--umiRange 24 31 + '
            f'-f {index_file} '
            f'--ref {ref} '
            f'{string1}'
            f'-o {self.sample} --od {self.outdir} --outputReadAssignment '
        )

        Assemble.run.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)


@utils.add_log
def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='R1 reads from match step', required=True)
        parser.add_argument('--fq2', help='R2 reads from match step', required=True)

    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)
    parser.add_argument('--speed_up', help='speed assemble for TCR/BCR seq data', action='store_true')       








