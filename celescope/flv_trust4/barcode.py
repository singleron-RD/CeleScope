from collections import Counter
from xopen import xopen
import pysam

from celescope.tools.__init__ import PATTERN_DICT
from celescope.tools import barcode as tools_barcode
from celescope.tools import utils

class Barcode(tools_barcode.Barcode):
    """
    ## Features

    - Demultiplex barcodes and UMIs.
    - Reverse complement the barcode to match RNA library barcodes.
    - Only reads with barcodes oberserbed in matched RNA library are kept.
    - If there are more than 80,000 reads for any barcodes, the reads are downsampled.


    ## Output

    - `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
    the read name is `{barcode}_{UMI}_{read ID}`.
    - `01.barcode/{sample}_1.fq(.gz)` Write barcode and umi to R1 read(can be directly used as input file of TRUST4).
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        # check flv
        for chemistry in self.chemistry_list:
            if chemistry != 'flv':
                raise Exception('chemistry should be `flv`!')

        self.match_num = 0 # read number match with flv_rna
        self.match_cbs = set() # barcode number match with flv_rna
        self.barcode_read_Counter = Counter()
        self.match_barcodes, _ = utils.get_barcode_from_match_dir(args.match_dir)

        self.fh_fq1 = xopen(self.out_fq1, 'w')


    @utils.add_log
    def run(self):
        """
        Extract barcode and UMI from R1. Filter reads with 
            - invalid polyT
            - low quality in barcode and UMI
            - invalid inlinker
            - invalid barcode
            
        for every sample
            get chemistry
            get linker_mismatch_dict and barcode_mismatch_dict
            for every read in read1
                filter
                write valid R2 read to file
        """

        for i in range(self.fq_number):

            chemistry = self.chemistry_list[i]
            lowNum = int(self.lowNum)
            lowQual = int(self.lowQual)
            if chemistry == 'scopeV1':
                lowNum = min(0, lowNum)
                lowQual = max(10, lowQual)
                Barcode.run.logger.info(f'scopeV1: lowNum={lowNum}, lowQual={lowQual} ')
            # get linker and whitelist
            bc_pattern = PATTERN_DICT[chemistry]
            if (bc_pattern):
                linker_file, whitelist_file = Barcode.get_scope_bc(chemistry)
            else:
                bc_pattern = self.pattern
                linker_file = self.linker
                whitelist_file = self.whitelist
            if not bc_pattern:
                raise Exception("invalid bc_pattern!")

            pattern_dict = Barcode.parse_pattern(bc_pattern)

            bool_T = True if 'T' in pattern_dict else False
            bool_L = True if 'L' in pattern_dict else False
            bool_whitelist = (whitelist_file is not None) and whitelist_file != "None"
            C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

            if bool_whitelist:
                barcode_set_list, barcode_mismatch_list = Barcode.parse_whitelist_file(whitelist_file,
                                                                               n_mismatch=1, n_repeat=len(pattern_dict['C']))
            if bool_L:
                linker_set_list, linker_mismatch_list = Barcode.parse_linker_file(linker_file)

            with pysam.FastxFile(self.fq1_list[i], persist=False) as fq1, \
                    pysam.FastxFile(self.fq2_list[i], persist=False) as fq2:
                for entry1, entry2 in zip(fq1, fq2):
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    self.total_num += 1

                    # polyT filter
                    if bool_T and (not self.allowNoPolyT):
                        polyT = Barcode.get_seq_str(seq1, pattern_dict['T'])
                        if polyT.count('T') < tools_barcode.MIN_T:
                            self.no_polyT_num += 1
                            if self.nopolyT:
                                Barcode.fh1_without_polyT.write(
                                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                                Barcode.fh2_without_polyT.write(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue

                    # lowQual filter
                    C_U_quals_ascii = Barcode.get_seq_str(
                        qual1, pattern_dict['C'] + pattern_dict['U'])
                    # C_U_quals_ord = [ord(q) - 33 for q in C_U_quals_ascii]
                    if lowQual > 0 and Barcode.low_qual(C_U_quals_ascii, lowQual, lowNum):
                        self.lowQual_num += 1
                        continue

                    # linker filter
                    if bool_L and (not self.allowNoLinker):
                        seq_str = Barcode.get_seq_str(seq1, pattern_dict['L'])
                        bool_valid, bool_corrected, _ = Barcode.check_seq_mismatch(
                            [seq_str], linker_set_list, linker_mismatch_list)
                        if not bool_valid:
                            self.no_linker_num += 1
                            if self.noLinker:
                                Barcode.fh1_without_linker.write(
                                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                                Barcode.fh2_without_linker.write(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue
                        elif bool_corrected:
                            self.linker_corrected_num += 1

                    # barcode filter
                    seq_list = Barcode.get_seq_list(seq1, pattern_dict, 'C')
                    if bool_whitelist:
                        bool_valid, bool_corrected, corrected_seq = Barcode.check_seq_mismatch(
                            seq_list, barcode_set_list, barcode_mismatch_list)

                        if not bool_valid:
                            self.no_barcode_num += 1
                            continue
                        elif bool_corrected:
                            self.barcode_corrected_num += 1
                        cb = corrected_seq
                    else:
                        cb = "".join(seq_list)
                    # flv barcode is reverse complement of the flv_rna barcode
                    cb = utils.reverse_complement(cb)
                    self.barcode_read_Counter.update(cb)

                    self.clean_num += 1
                    self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                    self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])

                    if cb in self.match_barcodes:
                        self.match_num += 1
                        self.match_cbs.add(cb)
                        if self.barcode_read_Counter[cb] <= 80000:
                            umi = Barcode.get_seq_str(seq1, pattern_dict['U'])
                            qual = 'F' * len(cb + umi)
                            self.fh_fq2.write(f'@{cb}_{umi}_{self.total_num}\n{seq2}\n+\n{qual2}\n')
                            self.fh_fq1.write(f'@{cb}_{umi}_{self.total_num}\n{cb}{umi}\n+\n{qual}\n')
                        
            self.run.logger.info(self.fq1_list[i] + ' finished.')

        self.close_files()
        self.add_step_metrics()

    @utils.add_log
    def add_step_metrics(self):
        super().add_step_metrics()

        self.add_metric(
            name='Valid Matched Reads',
            value=self.match_num,
            total=self.total_num,
            help_info='reads match with flv_rna cell barcodes'
        )

        self.add_metric(
            name='Matched Barcodes',
            value=len(self.match_cbs),
            help_info='barcodes match with flv_rna'
        )



@utils.add_log
def barcode(args):
    with Barcode(args, display_title='Demultiplexing') as runner:
        runner.run()


def get_opts_barcode(parser, sub_program=True):
    tools_barcode.get_opts_barcode(parser, sub_program)
    if sub_program:
        parser.add_argument('--match_dir', help='Matched scRNA-seq directory.', required=True)
