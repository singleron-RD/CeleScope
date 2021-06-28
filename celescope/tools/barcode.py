"""barcode step."""

import glob
import os
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations, product

import pysam
from xopen import xopen

import celescope.tools.utils as utils
from celescope.tools.__init__ import __PATTERN_DICT__
from celescope.tools.step import Step, s_common

MIN_T = 10


def seq_ranges(seq, pattern_dict):
    # get subseq with intervals in arr and concatenate
    return ''.join([seq[x[0]:x[1]]for x in pattern_dict])


def get_seq_list(seq, pattern_dict, abbr):
    # get subseq list
    return (seq[item[0]: item[1]] for item in pattern_dict[abbr])


@utils.add_log
def parse_pattern(pattern):
    pattern_dict = defaultdict(list)
    p = re.compile(r'([CLUNT])(\d+)')
    tmp = p.findall(pattern)
    if not tmp:
        parse_pattern.logger.error(f'Invalid pattern: {pattern}')
        sys.exit()
    start = 0
    for item in tmp:
        end = start + int(item[1])
        pattern_dict[item[0]].append([start, end])
        start = end
    return pattern_dict


def get_scope_bc(chemistry):
    import celescope
    root_path = os.path.dirname(celescope.__file__)
    if chemistry == 'scopeV1':
        return None, None
    linker_f = glob.glob(f'{root_path}/data/chemistry/{chemistry}/linker*')[0]
    whitelist_f = f'{root_path}/data/chemistry/{chemistry}/bclist'
    return linker_f, whitelist_f


def read_fastq(f):
    """
    Deprecated!
    Return tuples: (name, sequence, qualities).
    qualities is a string and it contains the unmodified, encoded qualities.
    """
    i = 3
    for i, line in enumerate(f):
        if i % 4 == 0:
            assert line.startswith('@'), ("Line {0} in FASTQ file is expected to start with '@', "
                                          "but found {1!r}".format(i + 1, line[:10]))
            name = line.strip()[1:]
        elif i % 4 == 1:
            sequence = line.strip()
        elif i % 4 == 2:
            line = line.strip()
            assert line.startswith('+'), ("Line {0} in FASTQ file is expected to start with '+', "
                                          "but found {1!r}".format(i + 1, line[:10]))
        elif i % 4 == 3:
            qualities = line.rstrip('\n\r')
            yield name, sequence, qualities
    if i % 4 != 3:
        raise Exception("FASTQ file ended prematurely")


def ord2chr(q, offset=33):
    return chr(int(q) + offset)


def qual_int(char, offset=33):
    return ord(char) - offset


def low_qual(quals, minQ, num):
    # print(ord('/')-33)           14
    return True if len([q for q in quals if qual_int(q) < minQ]) > num else False


def check_seq(seq_file, pattern_dict, seq_abbr):
    length = 0
    for item in pattern_dict[seq_abbr]:
        start = item[0]
        end = item[1]
        length += end - start
    with open(seq_file, 'r') as fh:
        for seq in fh:
            seq = seq.strip()
            if seq == '':
                continue
            if len(seq) != length:
                raise Exception(
                    f'length of L in pattern ({length}) do not equal to length in {seq_file} ({len(seq)}) !')


def get_mismatch(seq, n_mismatch=1, bases='ACGTN'):
    '''
    choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product

    Return:
    all mismatch <= n_mismatch dict. Key: mismatch_seq, value: orig_seq
    '''
    seq_set = set()
    mismatch_dict = {}
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in product(*seq_locs):
            seq_set.add(''.join(poss))
    for mismatch_seq in seq_set:
        mismatch_dict[mismatch_seq] = seq
    return mismatch_dict


@utils.add_log
def get_all_mismatch(seq_list, n_mismatch=1):
    '''
    Return:
    all mismatch <= n_mismatch dict. Key: mismatch_seq, value: orig_seq from seq_file
    '''
    mismatch_dict = {}
    correct_set = set()

    for seq in seq_list:
        seq = seq.strip()
        if seq == '':
            continue
        correct_set.add(seq)
        mismatch_dict.update(get_mismatch(seq, n_mismatch=n_mismatch))

    return correct_set, mismatch_dict


def check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list):
    '''
    Return bool_valid, bool_corrected
    '''
    bool_valid = True
    bool_corrected = False
    corrected_seq = ''
    for index, seq in enumerate(seq_list):
        if seq not in correct_set_list[index]:
            if seq not in mismatch_dict_list[index]:
                bool_valid = False
                return bool_valid, bool_corrected, corrected_seq
            else:
                bool_corrected = True
                corrected_seq += mismatch_dict_list[index][seq]
        else:
            corrected_seq += seq
    return bool_valid, bool_corrected, corrected_seq


class Chemistry():
    """
    Auto detect chemistry from read 1
    """

    def __init__(self, fq1):
        self.fq1 = fq1
        self.fq1_list = fq1.split(',')
        self.nRead = 10000

    @utils.add_log
    def check_chemistry(self):
        chemistry_list = []
        for fastq1 in self.fq1_list:
            print(fastq1)
            chemistry = self.get_chemistry(fastq1)
            chemistry_list.append(chemistry)
        if len(set(chemistry_list)) != 1:
            Chemistry.check_chemistry.logger.warning('multiple chemistry found!' + str(chemistry_list))
        return chemistry_list

    @utils.add_log
    def get_chemistry(self, fq1):
        '''
        'scopeV2.0.1': 'C8L16C8L16C8L1U8T18'
        'scopeV2.1.1': 'C8L16C8L16C8L1U12T18'
        'scopeV2.2.1': 'C8L16C8L16C8L1U12T18' with 4 types of linkers
        '''
        # init
        linker_4_file, _whitelist = get_scope_bc('scopeV2.2.1')
        linker_4_list, _num = utils.read_one_col(linker_4_file)
        linker_4_dict = defaultdict(int)
        linker_wrong_dict = defaultdict(int)
        pattern_dict = parse_pattern('C8L16C8L16C8L1U12T18')
        T4_n = 0
        L57C_n = 0

        with pysam.FastxFile(fq1) as fh:
            for _ in range(self.nRead):
                entry = fh.__next__()
                seq = entry.sequence
                L57C = seq[56]
                if L57C == 'C':
                    L57C_n += 1
                T4 = seq[65:69]
                if T4 == 'TTTT':
                    T4_n += 1
                linker = seq_ranges(seq, pattern_dict=pattern_dict['L'])
                if linker in linker_4_list:
                    linker_4_dict[linker] += 1
                else:
                    linker_wrong_dict[linker] += 1

        percent_T4 = T4_n / self.nRead
        percent_L57C = L57C_n / self.nRead
        Chemistry.get_chemistry.logger.info(f'percent T4: {percent_T4}')
        Chemistry.get_chemistry.logger.info(f'percent L57C: {percent_L57C}')
        if percent_T4 > 0.5:
            chemistry = 'scopeV2.0.1'
        else:
            # V2.1.1 or V2.2.1 or failed
            valid_linker_type = 0
            for linker in linker_4_dict:
                linker_4_dict[linker] = linker_4_dict[linker] / self.nRead
                if linker_4_dict[linker] > 0.05:
                    valid_linker_type += 1
            Chemistry.get_chemistry.logger.info(linker_4_dict)
            if valid_linker_type == 0:
                print(linker_wrong_dict)
                raise Exception(
                    'Auto chemistry detection failed! '
                    'If the sample is from Singleron, ask the technical staff you are connecting with for the chemistry used. '
                    'You need to use `--chemistry scopeV1` for scopeV1, and `--chemistry auto` should be fine for scopeV2.* '
                )
            elif valid_linker_type == 1:
                chemistry = 'scopeV2.1.1'
            elif valid_linker_type < 4:
                chemistry = 'scopeV2.2.1'
                Chemistry.get_chemistry.logger.warning(
                    f'chemistry scopeV2.2.1 only has {valid_linker_type} linker types!')
            else:
                chemistry = 'scopeV2.2.1'
        Chemistry.get_chemistry.logger.info(f'chemistry: {chemistry}')
        return chemistry


class Barcode(Step):
    """
    Features

    - Demultiplex barcodes.
    - Filter invalid R1 reads, which includes:
        - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
        - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
        - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
        - Low quality reads: low sequencing quality in barcode and UMI regions.

    Output

    - `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
    the read name is `{barcode}_{UMI}_{read ID}`.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        self.fq_number = len(self.fq1_list)
        if self.fq_number != len(self.fq2_list):
            raise Exception('fastq1 and fastq2 do not have same file number!')
        if args.chemistry == 'auto':
            ch = Chemistry(args.fq1)
            self.chemistry_list = ch.check_chemistry()
        else:
            self.chemistry_list = [args.chemistry] * self.fq_number
        self.barcode_corrected_num = 0
        self.linker_corrected_num = 0
        self.total_num = 0
        self.clean_num = 0
        self.no_polyT_num = 0
        self.lowQual_num = 0
        self.no_linker_num = 0
        self.no_barcode_num = 0
        self.barcode_qual_Counter = Counter()
        self.umi_qual_Counter = Counter()
        self.pattern = args.pattern
        self.linker = args.linker
        self.whitelist = args.whitelist
        self.lowNum = args.lowNum
        self.lowQual = args.lowQual
        self.allowNoPolyT = args.allowNoPolyT
        self.allowNoLinker = args.allowNoLinker
        self.nopolyT = args.nopolyT  # true == output nopolyT reads
        self.noLinker = args.noLinker

        # out file
        if args.gzip:
            suffix = ".gz"
        else:
            suffix = ""
        self.out_fq2 = f'{self.outdir}/{self.sample}_2.fq{suffix}'
        if self.nopolyT:
            self.nopolyT_1 = f'{self.outdir}/noPolyT_1.fq'
            self.nopolyT_2 = f'{self.outdir}/noPolyT_2.fq'
        if self.noLinker:
            self.noLinker_1 = f'{self.outdir}/noLinker_1.fq'
            self.noLinker_2 = f'{self.outdir}/noLinker_2.fq'

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

        fh3 = xopen(self.out_fq2, 'w')

        if self.nopolyT:
            fh1_without_polyT = xopen(self.nopolyT_1, 'w')
            fh2_without_polyT = xopen(self.nopolyT_2, 'w')

        if self.noLinker:
            fh1_without_linker = xopen(self.noLinker_1, 'w')
            fh2_without_linker = xopen(self.noLinker_2, 'w')

        for i in range(self.fq_number):

            chemistry = self.chemistry_list[i]
            lowNum = int(self.lowNum)
            Barcode.run.logger.info(f'lowQual score: {self.lowQual}')
            lowQual = int(self.lowQual)
            if chemistry == 'scopeV1':
                lowNum = min(0, lowNum)
                lowQual = max(10, lowQual)
                Barcode.run.logger.info(f'scopeV1: lowNum={lowNum}, lowQual={lowQual} ')
            # get linker and whitelist
            bc_pattern = __PATTERN_DICT__[chemistry]
            if (bc_pattern):
                (linker, whitelist) = get_scope_bc(chemistry)
            else:
                bc_pattern = self.pattern
                linker = self.linker
                whitelist = self.whitelist
            if not bc_pattern:
                raise Exception("invalid bc_pattern!")

            # parse pattern to dict, C8L10C8L10C8U8
            # defaultdict(<type 'list'>, {'C': [[0, 8], [18, 26], [36, 44]], 'U':
            # [[44, 52]], 'L': [[8, 18], [26, 36]]})
            pattern_dict = parse_pattern(bc_pattern)

            bool_T = True if 'T' in pattern_dict else False
            bool_L = True if 'L' in pattern_dict else False
            bool_whitelist = (whitelist is not None) and whitelist != "None"
            C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

            if bool_whitelist:
                seq_list, _ = utils.read_one_col(whitelist)
                barcode_correct_set, barcode_mismatch_dict = get_all_mismatch(seq_list, n_mismatch=1)
                barcode_correct_set_list = [barcode_correct_set] * 3
                barcode_mismatch_dict_list = [barcode_mismatch_dict] * 3
            if bool_L:
                seq_list, _ = utils.read_one_col(linker)
                check_seq(linker, pattern_dict, "L")
                linker_correct_set_list = []
                linker_mismatch_dict_list = []
                start = 0
                for item in pattern_dict['L']:
                    end = start + item[1] - item[0]
                    linker_seq_list = [seq[start:end] for seq in seq_list]
                    linker_correct_set, linker_mismatch_dict = get_all_mismatch(linker_seq_list, n_mismatch=2)
                    linker_correct_set_list.append(linker_correct_set)
                    linker_mismatch_dict_list.append(linker_mismatch_dict)
                    start = end

            fq1 = pysam.FastxFile(self.fq1_list[i], persist=False)
            fq2 = pysam.FastxFile(self.fq2_list[i], persist=False)

            for entry1 in fq1:
                entry2 = next(fq2)
                header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                self.total_num += 1

                # polyT filter
                if bool_T and (not self.allowNoPolyT):
                    polyT = seq_ranges(seq1, pattern_dict['T'])
                    if polyT.count('T') < MIN_T:
                        self.no_polyT_num += 1
                        if self.nopolyT:
                            fh1_without_polyT.write(
                                '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                            fh2_without_polyT.write(
                                '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                        continue

                # lowQual filter
                C_U_quals_ascii = seq_ranges(
                    qual1, pattern_dict['C'] + pattern_dict['U'])
                # C_U_quals_ord = [ord(q) - 33 for q in C_U_quals_ascii]
                if lowQual > 0 and low_qual(C_U_quals_ascii, lowQual, lowNum):
                    self.lowQual_num += 1
                    continue

                # linker filter
                if bool_L and (not self.allowNoLinker):
                    seq_list = get_seq_list(seq1, pattern_dict, 'L')
                    bool_valid, bool_corrected, _ = check_seq_mismatch(
                        seq_list, linker_correct_set_list, linker_mismatch_dict_list)
                    if not bool_valid:
                        self.no_linker_num += 1
                        if self.noLinker:
                            fh1_without_linker.write(
                                '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                            fh2_without_linker.write(
                                '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                        continue
                    elif bool_corrected:
                        self.linker_corrected_num += 1

                # barcode filter
                seq_list = get_seq_list(seq1, pattern_dict, 'C')
                if bool_whitelist:
                    bool_valid, bool_corrected, corrected_seq = check_seq_mismatch(
                        seq_list, barcode_correct_set_list, barcode_mismatch_dict_list)

                    if not bool_valid:
                        self.no_barcode_num += 1
                        continue
                    elif bool_corrected:
                        self.barcode_corrected_num += 1
                    cb = corrected_seq
                else:
                    cb = "".join(seq_list)

                umi = seq_ranges(seq1, pattern_dict['U'])

                self.clean_num += 1
                self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])

                fh3.write(f'@{cb}_{umi}_{self.total_num}\n{seq2}\n+\n{qual2}\n')
            Barcode.run.logger.info(self.fq1_list[i] + ' finished.')
        fh3.close()

        # logging
        Barcode.run.logger.info(
            f'processed reads: {utils.format_number(self.total_num)}. '
            f'valid reads: {utils.format_number(self.clean_num)}. '
        )

        Barcode.run.logger.info(f'no polyT reads number : {self.no_polyT_num}')
        Barcode.run.logger.info(f'low qual reads number: {self.lowQual_num}')
        Barcode.run.logger.info(f'no_linker: {self.no_linker_num}')
        Barcode.run.logger.info(f'no_barcode: {self.no_barcode_num}')
        Barcode.run.logger.info(f'corrected linker: {self.linker_corrected_num}')
        Barcode.run.logger.info(f'corrected barcode: {self.barcode_corrected_num}')

        if self.clean_num == 0:
            raise Exception(
                'no valid reads found! please check the --chemistry parameter.')

        # stat
        BarcodesQ30 = sum([self.barcode_qual_Counter[k] for k in self.barcode_qual_Counter if k >= ord2chr(
            30)]) / float(sum(self.barcode_qual_Counter.values())) * 100
        UMIsQ30 = sum([self.umi_qual_Counter[k] for k in self.umi_qual_Counter if k >= ord2chr(
            30)]) / float(sum(self.umi_qual_Counter.values())) * 100

        def cal_percent(x): return "{:.2%}".format((x + 0.0) / self.total_num)
        stat_info = '''
            Raw Reads: %s
            Valid Reads: %s(%s)
            Q30 of Barcodes: %.2f%%
            Q30 of UMIs: %.2f%%
        '''
        with open(self.stat_file, 'w') as fh:
            stat_info = stat_info % (utils.format_number(self.total_num), utils.format_number(self.clean_num),
                                     cal_percent(self.clean_num), BarcodesQ30,
                                     UMIsQ30)
            stat_info = re.sub(r'^\s+', r'', stat_info, flags=re.M)
            fh.write(stat_info)

        self.clean_up()


@utils.add_log
def barcode(args):
    step_name = "barcode"
    barcode_obj = Barcode(args, step_name)
    barcode_obj.run()


def get_opts_barcode(parser, sub_program=True):
    parser.add_argument(
        '--chemistry',
        help="""Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:  
- `auto` Default value. Used for Singleron GEXSCOPE libraries >= scopeV2 and automatically detects the combinations.  
- `scopeV1` Used for legacy Singleron GEXSCOPE scopeV1 libraries.  
- `customized` Used for user defined combinations. You need to provide `pattern`, `whitelist` and `linker` at the 
same time.""",
        choices=list(__PATTERN_DICT__.keys()),
        default='auto'
    )
    parser.add_argument(
        '--pattern',
        help="""The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  
- `C`: cell barcode  
- `L`: linker(common sequences)  
- `U`: UMI    
- `T`: poly T""",
    )
    parser.add_argument(
        '--whitelist',
        help='Cell barcode whitelist file path, one cell barcode per line.'
    )
    parser.add_argument(
        '--linker',
        help='Linker whitelist file path, one linker per line.'
    )
    parser.add_argument(
        '--lowQual',
        help='Default 0. Bases in cell barcode and UMI whose phred value are lower than \
lowQual will be regarded as low-quality bases.',
        type=int,
        default=0
    )
    parser.add_argument(
        '--lowNum',
        help='The maximum allowed lowQual bases in cell barcode and UMI.',
        type=int,
        default=2
    )
    parser.add_argument(
        '--nopolyT',
        help='Outputs R1 reads without polyT.',
        action='store_true',
    )
    parser.add_argument(
        '--noLinker',
        help='Outputs R1 reads without correct linker.',
        action='store_true',
    )
    parser.add_argument(
        '--allowNoPolyT',
        help="Allow valid reads without polyT.",
        action='store_true'
    )
    parser.add_argument(
        '--allowNoLinker',
        help="Allow valid reads without correct linker.",
        action='store_true'
    )
    parser.add_argument(
        '--gzip',
        help="Output gzipped fastq files.",
        action='store_true'
    )
    if sub_program:
        parser.add_argument('--fq1', help='R1 fastq file. Multiple files are separated by comma.', required=True)
        parser.add_argument('--fq2', help='R2 fastq file. Multiple files are separated by comma.', required=True)
        parser = s_common(parser)

    return parser
