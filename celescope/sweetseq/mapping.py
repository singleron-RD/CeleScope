import json

import pysam

from celescope.tools.count import Count
from celescope.tools import utils
from celescope.tools.barcode import Barcode
from celescope.tools.step import Step, s_common


def get_opts_mapping(parser, sub_program):
    parser.add_argument(
        "--fq_pattern",
        help="""R2 read pattern. The number after the letter represents the number of bases. The `fq_pattern` of CLindex is `L25C15`
`L` linker(common sequences)  
`C` tag barcode  
""",
        default='L25C15'
    )
    parser.add_argument(
        "--tag_barcode_fasta",
        help="""Required. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < len(tag barcode) / 10 + 1. 
If no such tag exists, the read is classified as invalid.
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
def mapping(args):
    with Mapping(args, display_title="Mapping") as runner:
        runner.run()


class Mapping(Step):
    """
    ## Features
    - Map R2 reads to the tag barcode.

    ## Output
    - raw_read_count.json: TODO
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # read args
        self.fq = args.fq
        self.fq_pattern = args.fq_pattern
        self.linker_fasta = args.linker_fasta
        self.tag_barcode_fasta = args.tag_barcode_fasta

        # process
        self.barcode_dict, self.barcode_length = utils.read_fasta(self.tag_barcode_fasta, equal=True)
        if len(self.barcode_dict) > 1:
            raise ValueError("More than one barcode in tag_barcode_fasta.")
        if self.linker_fasta and self.linker_fasta != 'None':
            self.linker_dict, self.linker_length = utils.read_fasta(self.linker_fasta, equal=True)
        else:
            self.linker_dict, self.linker_length = {}, 0
        self.pattern_dict = Barcode.parse_pattern(self.fq_pattern)

        # check barcode length
        barcode_start, barcode_end = self.pattern_dict["C"][0]
        # end - start
        pattern_barcode_length = barcode_end - barcode_start
        if pattern_barcode_length != self.barcode_length:
            raise Exception(
                f'''barcode fasta length {self.barcode_length} 
                != pattern barcode length {pattern_barcode_length}'''
            )

        self.read_count_dict = utils.genDict(dim=3)
        self.res_sum_dic = utils.genDict(dim=2)
        self.match_barcode = []

        # metrics
        self.total_reads = 0
        self.reads_unmapped_too_short = 0
        self.reads_unmapped_invalid_iinker = 0
        self.reads_unmapped_invalid_barcode = 0
        self.reads_mapped = 0
        self.raw_umi = 0
        self.total_corrected_umi = 0

        # out files
        self.raw_read_count_file = f'{self.out_prefix}_raw_read_count.json'
        self.corrected_read_count_file = f'{self.out_prefix}_corrected_read_count.json'
        self.UMI_count_file = f'{self.out_prefix}_UMI_count.csv'


    def map_read(self):

        with pysam.FastxFile(self.fq) as infile:
            for record in infile:
                self.total_reads += 1
                attr = str(record.name).strip("@").split("_")
                barcode = str(attr[0])
                umi = str(attr[1])
                seq = record.sequence

                if self.linker_length != 0:
                    seq_linker = Barcode.get_seq_str(seq, self.pattern_dict['L'])
                    if len(seq_linker) < self.linker_length:
                        self.reads_unmapped_too_short += 1
                        continue
                if self.barcode_dict:
                    seq_barcode = Barcode.get_seq_str(seq, self.pattern_dict['C'])
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
                    self.reads_unmapped_invalid_iinker += 1
                    continue

                # check barcode
                valid_barcode = False
                for barcode_name in self.barcode_dict:
                    if utils.hamming_correct(self.barcode_dict[barcode_name], seq_barcode):
                        self.read_count_dict[barcode][barcode_name][umi] += 1
                        valid_barcode = True
                        break

                if not valid_barcode:
                    self.reads_unmapped_invalid_barcode += 1
                    continue

                # mapped
                self.reads_mapped += 1

    def add_map_metrics(self):
    
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
            value=self.reads_unmapped_invalid_iinker,
            total=self.total_reads,
            help_info="Unmapped R2 reads because of too many mismatches in linker sequence"
        )
        self.add_metric(
            name='Reads Unmapped Invalid Barcode',
            value=self.reads_unmapped_invalid_barcode,
            total=self.total_reads,
            help_info="Unmapped R2 reads because of too many mismatches in tag-barcode sequence"
        )


    def write_raw_read_count_file(self):
        with open(self.raw_read_count_file, 'w') as fp:
            json.dump(self.read_count_dict, fp, indent=4)

    def correct_umi(self):
        for barcode in self.count_dict:
            for ref in self.count_dict[barcode]:
                self.raw_umi += len(self.count_dict[barcode][ref])
                n_corrected_umi, _n_corrected_read = Count.correct_umi(self.count_dict[barcode][ref])
                if self.debug:
                    print(f'{barcode} {ref} {n_corrected_umi}')
                self.total_corrected_umi += n_corrected_umi

    def add_correct_umi_metrics(self):
        self.add_metric(
            name='Number of Raw UMI',
            value=self.raw_umi,
            help_info='number of total raw UMI',
        )

        self.add_metric(
            name='Number of Corrected UMI',
            value=self.total_corrected_umi,
            total=self.raw_umi,
            help_info='correct sequencing errors in the UMI sequences ',
        )

    def write_corrected_read_count_file(self):
        with open(self.corrected_read_count_file, 'w') as fp:
            json.dump(self.read_count_dict, fp, indent=4)

    def write_UMI_count_file(self):
        pass

        
    @utils.add_log
    def run(self):
        self.map_read()
        self.write_raw_read_count_file()
        self.add_map_metrics()
        '''
        if not self.args.not_correct_umi:
            self.correct_umi()
            self.write_corrected_read_count_file()
            self.add_correct_umi_metrics()
        self.write_tsne_UMI_file()
        '''
