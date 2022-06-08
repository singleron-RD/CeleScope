from celescope.tools import utils
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


@utils.add_log
def summarize(args):
    with Summarize(args, display_title="Cells") as runner:
        runner.run()


def get_opts_summarize(parser, sub_program):
    tools_summarize.get_opts_summarize(parser, sub_program)

