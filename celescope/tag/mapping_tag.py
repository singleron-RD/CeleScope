"""
map read2 to barcode_fasta
"""

import pandas as pd
import pysam

from celescope.tools import utils
from celescope.tools.barcode import Barcode
from celescope.tools.step import Step, s_common

# n_mismatch = 1 if n_tag_barcode > N_TAG_BARCODE_THRESHOLD else 2
N_TAG_BARCODE_THRESHOLD = 10000

def get_opts_mapping_tag(parser, sub_program):
    parser.add_argument(
        "--fq_pattern",
        help="""R2 read pattern. The number after the letter represents the number of bases. The `fq_pattern` of CLindex is `L25C15`
`L` linker(common sequences)  
`C` tag barcode  
""",
        default='L25C15'
    )
    parser.add_argument(
        "--barcode_fasta",
        help="""Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < 2. 
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
        self.pattern_dict = Barcode.parse_pattern(self.fq_pattern)

        self.barcode_dict, self.barcode_length = utils.read_fasta(self.barcode_fasta, equal=True)
        len_C = Barcode.get_abbr_len(self.pattern_dict, 'C')
        if len_C != self.barcode_length:
            raise ValueError(f"""The length of tag barcode in fq_pattern({len_C}) != 
                length of tag barcode in barcode_fasta({self.barcode_length})""")

        if self.linker_fasta and self.linker_fasta != 'None':
            self.linker_dict, self.linker_length = utils.read_fasta(self.linker_fasta, equal=True)
            len_L = Barcode.get_abbr_len(self.pattern_dict, 'L')
            if len_L != self.linker_length:
                raise ValueError(f"""The length of linker in fq_pattern({len_L}) != 
                    length of linker in linker_fasta({self.linker_length})""")
        else:
            self.linker_dict, self.linker_length = {}, 0




        # check barcode length
        barcode1 = self.pattern_dict["C"][0]
        # end - start
        pattern_barcode_length = barcode1[1] - barcode1[0]
        if pattern_barcode_length != self.barcode_length:
            raise Exception(
                f'''barcode fasta length {self.barcode_length} 
                != pattern barcode length {pattern_barcode_length}'''
            )

        # mismatch
        self.mismatch_dict = self.get_tag_barcode_mismatch_dict()

        # variables
        self.total_reads = 0
        self.reads_unmapped_too_short = 0
        self.reads_unmapped_invalid_linker = 0
        self.reads_unmapped_invalid_barcode = 0
        self.reads_mapped = 0
        self.res_dic = utils.genDict()
        self.res_sum_dic = utils.genDict(dim=2)
        self.match_barcode = []
        self.invalid_barcode_dict = utils.genDict(dim=1)

        # out files
        self.read_count_file = f'{self.outdir}/{self.sample}_read_count.tsv'
        self.UMI_count_file = f'{self.outdir}/{self.sample}_UMI_count.tsv'
        self.invalid_barcode_file = f'{self.outdir}/{self.sample}_invalid_barcode.tsv'
        self.stat_file = f'{self.outdir}/stat.txt'

    @utils.add_log
    def get_tag_barcode_mismatch_dict(self):
        mismatch_dict = {}
        n_mismatch = 1 if len(self.barcode_dict) > N_TAG_BARCODE_THRESHOLD else 2
        for seq_id, seq in self.barcode_dict.items():
            for mismatch_seq in Barcode.findall_mismatch(seq, n_mismatch=n_mismatch):
                mismatch_dict[mismatch_seq] = seq_id   
        
        return mismatch_dict


    def check_barcode_with_mismatch(self, barcode, seq_barcode, umi):
        """
        Args:
            barcode: cell barcode
            seq_barcode: tag barcode sequence
            umi: UMI sequence
        """
        if seq_barcode in self.mismatch_dict:
            seq_id = self.mismatch_dict[seq_barcode]
            self.res_dic[barcode][seq_id][umi] += 1
            self.reads_mapped += 1
        else:
            self.reads_unmapped_invalid_barcode += 1
            self.invalid_barcode_dict[seq_barcode] += 1

    def process_read(self):

        with pysam.FastxFile(self.fq) as infile:
            for record in infile:
                self.total_reads += 1
                attr = str(record.name).strip("@").split("_")
                barcode = str(attr[0])
                umi = str(attr[1])
                seq = record.sequence

                if self.linker_length != 0:
                    seq_linker = Barcode.get_seq_str_no_exception(seq, self.pattern_dict['L'])
                    if len(seq_linker) < self.linker_length:
                        self.reads_unmapped_too_short += 1
                        continue
                if self.barcode_dict:
                    seq_barcode = Barcode.get_seq_str_no_exception(seq, self.pattern_dict['C'])
                    if self.barcode_length != len(seq_barcode):
                        miss_length = self.barcode_length - len(seq_barcode)
                        if miss_length > 2:
                            self.reads_unmapped_too_short += 1
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
                    self.reads_unmapped_invalid_linker += 1
                    continue

                # check barcode
                self.check_barcode_with_mismatch(barcode, seq_barcode, umi)

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

        # write invalid seq_barcode to file
        with open(self.invalid_barcode_file, "w") as f:
            for seq_barcode,count in sorted(self.invalid_barcode_dict.items(), key=lambda x: x[1], reverse=True):
                f.write(f"{seq_barcode}\t{count}\n")

        # add metrics
        self.add_metric(
            name='Reads Mapped',
            value=self.reads_mapped,
            total=self.total_reads,
            help_info="R2 reads that successfully mapped to linker and tag-barcode"
        )
        self.add_metric(
            name='Reads Unmapped too Short',
            value=self.reads_unmapped_too_short,
            total=self.total_reads,
            help_info="Unmapped R2 reads because read length < linker length + tag-barcode length"
        )
        self.add_metric(
            name='Reads Unmapped Invalid Linker',
            value=self.reads_unmapped_invalid_linker,
            total=self.total_reads,
            help_info="Unmapped R2 reads because of too many mismatches in linker sequence"
        )
        self.add_metric(
            name='Reads Unmapped Invalid Barcode',
            value=self.reads_unmapped_invalid_barcode,
            total=self.total_reads,
            help_info="Unmapped R2 reads because of too many mismatches in tag-barcode sequence"
        )

    @utils.add_log
    def run(self):
        self.process_read()
