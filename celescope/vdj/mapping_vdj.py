'''
vdj mapping
'''

from tkinter.tix import InputOnly
import pandas as pd
import pysam

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj.__init__ import CHAINS
from celescope.__init__ import HELP_DICT


class Mapping_vdj(Step):
    """
    ## Features
    - Align R2 reads to IGMT(http://www.imgt.org/) database sequences with mixcr.

    ## Output
    - `{sample}_consensus.fasta` Fasta file after UMI consensus.

    - `{sample}_UMI_count_unfiltered.tsv` UMI reading for each (barcode, chain, VJ_pair) combination.

    - `{sample}_UMI_count_filtered.tsv` For each (barcode, chain) combination, only the record with the 
    most VJ_pair UMI reads is kept.

    - `{sample}_align.txt` Result report.

    - `{sample}_alignments.txt` The alignment result of each UMI/read.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # set
        self.read_type = "UMIs"
        if args.not_consensus:
            self.read_type = 'Reads'
        self.chains = CHAINS[args.type]

        # out files
        self.UMI_count_unfiltered_file = f'{self.out_prefix}_UMI_count_unfiltered.tsv'
        self.UMI_count_filtered_file = f'{self.out_prefix}_UMI_count_filtered.tsv'
        self.mixcr_report = f"{self.out_prefix}_align.txt"
        self.not_align_fq = f"{self.out_prefix}_not_align.fq"
        self.read2_vdjca = f"{self.out_prefix}_read2.vdjca"
        self.alignments = f"{self.out_prefix}_alignments.txt"

    @utils.add_log
    def run_mixcr(self):
        cmd = (
            'mixcr align '
            '--force-overwrite '
            f'--species {self.args.species} '
            f'-t {self.args.thread} '
            f'--not-aligned-R1 {self.not_align_fq} '
            f'--report {self.mixcr_report} '
            '-OallowPartialAlignments=true '
            '-OvParameters.geneFeatureToAlign=VTranscriptWithP '
            f'{self.args.fq} {self.read2_vdjca} '
            '\n'
            'mixcr exportAlignments '
            f'{self.read2_vdjca} {self.alignments} '
            '-readIds --force-overwrite -vGene -dGene -jGene -cGene '
            '-nFeature CDR3 -aaFeature CDR3 '
        )

        self.debug_subprocess_call(cmd)

    @utils.add_log
    def mixcr_summary(self, total_read, df_align):

        self.add_help_content(
            name='',
            content='If `--not_consensus` argument was used, reads were used instead of UMIs.',
        )

        align_read = df_align.shape[0]
        self.add_metric(
            name=f"{self.read_type} Mapped to Any VDJ Gene",
            value=align_read,
            total=total_read,
            help_info=f"{self.read_type} Mapped to any germline VDJ gene segments"
        )

        # CDR3
        df_CDR3 = df_align[~pd.isnull(df_align["aaSeqCDR3"])]
        align_read_with_CDR3 = df_CDR3.shape[0]
        self.add_metric(
            name=f"{self.read_type} with CDR3",
            value=align_read_with_CDR3,
            total=total_read,
            help_info=f"{self.read_type} with CDR3 sequence"
        )

        # correct CDR3
        df_correct_CDR3 = df_CDR3[~(df_CDR3["aaSeqCDR3"].str.contains(r"\*"))]
        align_read_with_correct_CDR3 = df_correct_CDR3.shape[0]
        self.add_metric(
            name=f"{self.read_type} with Correct CDR3",
            value=align_read_with_correct_CDR3,
            total=total_read,
            help_info=f"{self.read_type} with CDR3 might have stop codon and these UMIs(or Reads) are classified as incorrect"
        )

        # VDJ
        df_VJ = df_correct_CDR3[
            (~pd.isnull(df_correct_CDR3['bestVGene'])) &
            (~pd.isnull(df_correct_CDR3['bestJGene']))
        ]
        df_VJ = df_VJ[df_VJ.bestVGene.str[:3] == df_VJ.bestJGene.str[:3]]
        df_VJ["chain"] = df_VJ.bestVGene.str[:3]
        df_VJ["VJ_pair"] = df_VJ["bestVGene"] + "_" + df_VJ["bestJGene"]
        Reads_Mapped_Confidently_to_VJ_Gene = df_VJ.shape[0]
        self.add_metric(
            name=f"{self.read_type} Mapped Confidently to VJ Gene",
            value=Reads_Mapped_Confidently_to_VJ_Gene,
            total=total_read,
            help_info=f"{self.read_type} mapped to VJ gene pairs and with correct CDR3"
        )

        # chain
        for chain in self.chains:
            df_chain = df_VJ[df_VJ.chain == chain]
            Reads_Mapped_to_chain = df_chain.shape[0]
            self.add_metric(
                name=f"{self.read_type} Mapped to {chain}",
                value=Reads_Mapped_to_chain,
                total=total_read,
                help_info=f"{self.read_type} mapped confidently to {chain}"
            )

        # unique UMI
        df_UMI = df_VJ.drop_duplicates(subset=["barcode", "UMI"], keep="first")

        # filter1: keep top 1 in each combinations
        groupby_elements = [
            'barcode',
            'chain',
            'bestVGene',
            'bestJGene',
            'aaSeqCDR3',
            'nSeqCDR3',
        ]
        df_UMI_count = df_UMI.groupby(
            groupby_elements, as_index=False).agg({"UMI": "count"})
        df_UMI_count = df_UMI_count.sort_values("UMI", ascending=False)
        # out unfiltered
        df_UMI_count.to_csv(self.UMI_count_unfiltered_file, sep="\t", index=False)

        df_UMI_count_filter = df_UMI_count.groupby(
            ["barcode", "chain"], as_index=False).head(1)
        # out filtered
        df_UMI_count_filter.to_csv(
            self.UMI_count_filtered_file,
            sep="\t",
            index=False
        )

        if self.args.debug:
            unique_UMI = df_UMI.shape[0]
            self.add_metric(
                name="UMI unique count",
                value=unique_UMI,
                total=align_read_with_correct_CDR3,
            )

            UMI_after_Contamination_Filtering = df_UMI_count.filter.UMI.sum()
            self.add_metric(
                name="UMI after Contamination Filtering",
                value=UMI_after_Contamination_Filtering,
                total=unique_UMI,
            )

    @utils.add_log
    def fastq_to_dataframe(self):
        # read input_file
        with pysam.FastxFile(self.args.fq) as fh:
            index = 0
            read_row_list = []
            for entry in fh:
                attr = entry.name.split("_")
                barcode = attr[0]
                umi = attr[1]
                dic = {"readId": index, "barcode": barcode, "UMI": umi}
                read_row_list.append(dic)
                index += 1
            df_fastq = pd.DataFrame(read_row_list, columns=["readId", "barcode", "UMI"])
        return df_fastq

    def get_df_align(self, df_fastq):
        alignments = pd.read_csv(self.alignments, sep="\t")
        alignments.readId = alignments.readId.astype(int)
        df_fastq.readId = df_fastq.readId.astype(int)
        df_align = pd.merge(df_fastq, alignments, on="readId", how="right")
        return df_align

    def run(self):
        self.run_mixcr()
        df_fastq = self.fastq_to_dataframe()
        total_read = df_fastq.shape[0]
        df_align = self.get_df_align(df_fastq)
        self.mixcr_summary(total_read, df_align)
        self._clean_up()


@utils.add_log
def mapping_vdj(args):
    # TODO
    # add TCR or BCR prefix to distinguish them in html report summary; should improve
    with Mapping_vdj(args, display_title="Mapping") as runner:
        runner.run()


def get_opts_mapping_vdj(parser, sub_program):
    parser.add_argument("--type", help=HELP_DICT['type'], required=True)
    parser.add_argument(
        '--species',
        choices=['hs', 'mmu'],
        help=HELP_DICT['species'],
        default='hs'
    )
    parser.add_argument("--not_consensus", action='store_true', help=HELP_DICT['not_consensus'])
    if sub_program:
        parser.add_argument("--fq",help=HELP_DICT['fq'],required=True,)
        parser = s_common(parser)
