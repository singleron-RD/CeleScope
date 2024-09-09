from celescope.tools import consensus as supert_cs
from celescope.tools import utils
from xopen import xopen
import pysam


class Consensus(supert_cs.Consensus):
    """
    ## Features
    - Consensus all the reads of the same (barcode, UMI) combinations into one read(UMI). It will go through the sequence residue by residue and
    count up the number of each type of residue (ie. A or G or T or C for DNA) in all sequences in the
    alignment. If the following conditions are met, the consensus sequence will be the most common residue in the alignment:
    1. the percentage of the most common residue type > threshold(default: 0.5);
    2. most common residue reads >= min_consensus_read;
    otherwise an ambiguous character(N) will be added.

    ## Output
    - `{sample}_consensus.fq` Fastq file after consensus.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.filtered_consensus_out = f"{self.out_prefix}_filtered_consensus.fasta"

    def __exit__(self, *args, **kwargs):
        self.filter_consensus_out()
        self._clean_up()

    @utils.add_log
    def filter_consensus_out(self):
        out_h = xopen(self.filtered_consensus_out, "w")
        n_umi_count = 0
        n_filter_umi_count = 0

        with pysam.FastxFile(self.consensus_out) as f:
            for read in f:
                n_umi_count += 1
                # if set(read.sequence) == {'N'}:
                # if read.sequence.count('N') > 10:
                if read.sequence.count("N") > 5:
                    n_filter_umi_count += 1
                else:
                    out_h.write(utils.fasta_line(read.name, read.sequence))
        out_h.close()

        self.add_metric(
            name="Filtered UMI Counts",
            value=n_umi_count - n_filter_umi_count,
            total=n_umi_count,
            help_info="Filtered UMI from Consensus UMI fasta",
        )


@utils.add_log
def consensus(args):
    if args.not_consensus:
        return
    with Consensus(args, display_title="Consensus") as runner:
        runner.run()


def get_opts_consensus(parser, sub_program):
    supert_cs.get_opts_consensus(parser, sub_program)
