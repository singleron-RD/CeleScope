"""
map read2 to barcode_fasta
"""

import pandas as pd
import pysam

from celescope.tools import utils
import celescope.tools.barcode as Barcode
from celescope.tools.barcode import parse_pattern
from celescope.tools.step import Step, s_common


def get_opts_mapping_tag(parser, sub_program):
    parser.add_argument(
        "--fq_pattern",
        help="""Required. R2 read pattern. The number after the letter represents the number of bases.         
`L` linker(common sequences)  
`C` tag barcode  
""",
        required=True
    )
    parser.add_argument(
        "--barcode_fasta",
        help="""Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < len(tag barcode) / 10 + 1. 
If no such tag exists, the read is classified as invalid.

You can find the barcode fasta file under `celescope/data/Clindex`
```
>CLindex_TAG_1
CGTGTTAGGGCCGAT
>CLindex_TAG_2
GAGTGGTTGCGCCAT
>CLindex_TAG_3
AAGTTGCCAAGGGCC
>CLindex_TAG_4
TAAGAGCCCGGCAAG
>CLindex_TAG_5
TGACCTGCTTCACGC
>CLindex_TAG_6
GAGACCCGTGGAATC
>CLindex_TAG_7
GTTATGCGACCGCGA
>CLindex_TAG_8
ATACGCAGGGTCCGA
>CLindex_TAG_9
AGCGGCATTTGGGAC
>CLindex_TAG_10
TCGCCAGCCAAGTCT
>CLindex_TAG_11
ACCAATGGCGCATGG
>CLindex_TAG_12
TCCTCCTAGCAACCC
>CLindex_TAG_13
GGCCGATACTTCAGC
>CLindex_TAG_14
CCGTTCGACTTGGTG
>CLindex_TAG_15
CGCAAGACACTCCAC
>CLindex_TAG_16
CTGCAACAAGGTCGC
```
""",
        required=True,
    )
    parser.add_argument(
        "--linker_fasta",
        help="""Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.
""",
    )
    if sub_program:
        s_common(parser)
        parser.add_argument("--fq", help="R2 read fastq.", required=True)


@utils.add_log
def mapping_tag(args):
    with Mapping_tag(args, display_title="Mapping") as runner:
        runner.run()


class Mapping_tag(Step):
    """
    ## Features
    - Align R2 reads to the tag barcode fasta.

    ## Output

    - `{sample}_read_count.tsv` tab-delimited text file with 4 columns.

        `barcode` cell barcode  
        `tag_name`  tag name in barcode_fasta  
        `UMI`   UMI sequence  
        `read_count` read count per UMI  
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # read args
        self.fq = args.fq
        self.fq_pattern = args.fq_pattern
        self.linker_fasta = args.linker_fasta
        self.barcode_fasta = args.barcode_fasta

        # process
        self.barcode_dict, self.barcode_length = utils.read_fasta(self.barcode_fasta, equal=True)
        if self.linker_fasta and self.linker_fasta != 'None':
            self.linker_dict, self.linker_length = utils.read_fasta(self.linker_fasta, equal=True)
        else:
            self.linker_dict, self.linker_length = {}, 0
        self.pattern_dict = parse_pattern(self.fq_pattern)

        # check barcode length
        barcode1 = self.pattern_dict["C"][0]
        # end - start
        pattern_barcode_length = barcode1[1] - barcode1[0]
        if pattern_barcode_length != self.barcode_length:
            raise Exception(
                f'''barcode fasta length {self.barcode_length} 
                != pattern barcode length {pattern_barcode_length}'''
            )

        self.res_dic = utils.genDict()
        self.res_sum_dic = utils.genDict(dim=2)
        self.match_barcode = []

        # out files
        self.read_count_file = f'{self.outdir}/{self.sample}_read_count.tsv'
        self.UMI_count_file = f'{self.outdir}/{self.sample}_UMI_count.tsv'
        self.stat_file = f'{self.outdir}/stat.txt'

    def process_read(self):
        total_reads = 0
        reads_unmapped_too_short = 0
        reads_unmapped_invalid_iinker = 0
        reads_unmapped_invalid_barcode = 0
        reads_mapped = 0

        with pysam.FastxFile(self.fq) as infile:
            for record in infile:
                total_reads += 1
                attr = str(record.name).strip("@").split("_")
                barcode = str(attr[0])
                umi = str(attr[1])
                seq = record.sequence

                if self.linker_length != 0:
                    seq_linker = Barcode.get_seq_str(seq, self.pattern_dict['L'])
                    if len(seq_linker) < self.linker_length:
                        reads_unmapped_too_short += 1
                        continue
                if self.barcode_dict:
                    seq_barcode = Barcode.get_seq_str(seq, self.pattern_dict['C'])
                    if self.barcode_length != len(seq_barcode):
                        miss_length = self.barcode_length - len(seq_barcode)
                        if miss_length > 2:
                            reads_unmapped_too_short += 1
                            continue
                        seq_barcode = seq_barcode + "A" * miss_length

                # check linker
                if self.linker_length != 0:
                    valid_linker = False
                    for linker_name in self.linker_dict:
                        if utils.hamming_correct(self.linker_dict[linker_name], seq_linker):
                            valid_linker = True
                            break
                else:
                    valid_linker = True

                if not valid_linker:
                    reads_unmapped_invalid_iinker += 1
                    continue

                # check barcode
                valid_barcode = False
                for barcode_name in self.barcode_dict:
                    if utils.hamming_correct(self.barcode_dict[barcode_name], seq_barcode):
                        self.res_dic[barcode][barcode_name][umi] += 1
                        valid_barcode = True
                        break

                if not valid_barcode:
                    reads_unmapped_invalid_barcode += 1
                    continue

                # mapped
                reads_mapped += 1

        # write dic to pandas df
        rows = []
        for barcode in self.res_dic:
            for tag_name in self.res_dic[barcode]:
                for umi in self.res_dic[barcode][tag_name]:
                    rows.append([barcode, tag_name, umi,
                                 self.res_dic[barcode][tag_name][umi]])
        df_read_count = pd.DataFrame(rows)
        df_read_count.rename(
            columns={
                0: "barcode",
                1: "tag_name",
                2: "UMI",
                3: "read_count"
            }, inplace=True)
        df_read_count.to_csv(
            self.read_count_file, sep="\t", index=False)

        # add metrics
        self.add_metric(
            name='Reads Mapped',
            value=reads_mapped,
            total=total_reads,
            help_info="R2 reads that successfully mapped to linker and tag-barcode"
        )
        self.add_metric(
            name='Reads Unmapped too Short',
            value=reads_unmapped_too_short,
            total=total_reads,
            help_info="Unmapped R2 reads because read length < linker length + tag-barcode length"
        )
        self.add_metric(
            name='Reads Unmapped Invalid Linker',
            value=reads_unmapped_invalid_iinker,
            total=total_reads,
            help_info="Unmapped R2 reads because of too many mismatches in linker sequence"
        )
        self.add_metric(
            name='Reads Unmapped Invalid Barcode',
            value=reads_unmapped_invalid_barcode,
            total=total_reads,
            help_info="Unmapped R2 reads because of too many mismatches in tag-barcode sequence"
        )

    @utils.add_log
    def run(self):
        self.process_read()
