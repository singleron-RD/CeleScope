from celescope.tools import utils
from celescope.flv_trust4.summarize import *
import celescope.flv_trust4.summarize as tools_summarize

class Summarize(tools_summarize.Summarize):

    """
    ## Features

    - TCR/BCR full length assembly results.

    ## Output
    - `04.summarize/clonetypes.csv` High-level descriptions of each clonotype.
    - `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
    - `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
    - `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig.csv.
    - `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
    """
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.seqtype = args.seqtype
        self.fq2 = args.fq2
        self.diffuseFrac = args.diffuseFrac
        self.assembled_fa = f'{args.assemble_out}/{self.sample}_assembled_reads.fa'
        self.trust_report = f'{args.assemble_out}/reportfl.tsv'
        self.annot = f'{args.assemble_out}/{self.sample}_annot.fa'

        # if --diffuseFrac provided
        if self.diffuseFrac:
            self.barcode_report = f'{args.assemble_out}/barcoderepfl.tsv'
        else:
            self.barcode_report = f'{args.assemble_out}/barcoderep.tsv'
        
        self.coef = int(args.coef)
        self.target_weight = args.target_weight
        if args.target_cell_barcode:
            self.target_barcodes, self.expected_target_cell_num = utils.read_one_col(args.target_cell_barcode)
        else:
            self.target_barcodes = None
            self.expected_target_cell_num = args.expected_target_cell_num

        self.matrix_file = utils.get_matrix_dir_from_match_dir(args.match_dir)
        self.chains, self.paired_groups = self._parse_seqtype(self.seqtype)
        self.record_file = f'{self.outdir}/Cell_num.txt'


@utils.add_log
def summarize(args):
    with Summarize(args, display_title="Cells") as runner:
        runner.run()


def get_opts_summarize(parser, sub_program):
    tools_summarize.get_opts_summarize(parser, sub_program)

