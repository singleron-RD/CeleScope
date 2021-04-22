import os
import re
import io
import gzip
import subprocess
import sys
import glob
import pandas as pd
import pysam

from collections import defaultdict, Counter
from itertools import combinations, permutations, islice
from xopen import xopen

import celescope.tools.utils as utils
from celescope.tools.__init__ import __PATTERN_DICT__
from celescope.tools.Chemistry import Chemistry
from celescope.tools.Step import Step, s_common


def seq_ranges(seq, pattern_dict):
    # get subseq with intervals in arr and concatenate
    return ''.join([seq[x[0]:x[1]]for x in pattern_dict])


def generate_mis_seq(seq, n=1, bases='ACGTN'):
    # 以随机bases中的碱基替换seq中的n个位置，产生的错配字典
    # 返回字典，错配序列为key，
    # (正确序列，错配碱基数目，错配碱基位置，原始碱基，新碱基)组成的元组
    # 作为字典的值

    length = len(seq)
    assert length >= n, "err number should not be larger than sequence length!"
    res = {}
    seq_arr = list(seq)
    pos_group = list(combinations(range(0, length), n))
    bases_group = list(permutations(bases, n))

    for g in pos_group:
        for b in bases_group:
            seq_tmp = seq_arr[:]
            mis_num = n
            raw_tmp = []
            for i in range(n):
                raw_base = seq_tmp[g[i]]
                new_base = b[i]

                if raw_base == new_base:
                    mis_num -= 1

                raw_tmp.append(raw_base)
                seq_tmp[g[i]] = new_base

            if mis_num != 0:
                res[''.join(seq_tmp)] = (seq, mis_num, ','.join(
                    [str(i) for i in g]), ','.join(raw_tmp), ','.join(b))
    return(res)


@utils.add_log
def generate_seq_dict(seqlist, n=1):
    seq_dict = {}
    with open(seqlist, 'r') as fh:
        for seq in fh:
            seq = seq.strip()
            if seq == '':
                continue
            seq_dict[seq] = (seq, 0, -1, 'X', 'X')
            for k, v in generate_mis_seq(seq, n).items():
                # duplicate key
                if k in seq_dict:
                    generate_seq_dict.logger.warning('barcode %s, %s\n%s, %s' %
                                                     (v, k, seq_dict[k], k))
                else:
                    seq_dict[k] = v
    return seq_dict


@utils.add_log
def parse_pattern(pattern):
    # 解析接头结构，返回接头结构字典
    # key: 字母表示的接头, value: 碱基区间列表
    # eg.: C8L10C8L10C8U8T30
    # defaultdict(<type 'list'>:
    # {'C': [[0, 8], [18, 26], [36, 44]], 'U': [[44, 52]], 'L': [[8, 18], [26, 36]], 'T': [[52, 82]]})
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


def no_polyT(seq, strictT=0, minT=10):
    # strictT requires the first nth position to be T
    if seq[:strictT] != 'T' * strictT or seq.count('T') < minT:
        return True
    else:
        return False


def no_linker(seq, linker_dict):
    return False if seq in linker_dict else True


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


class Barcode(Step):

    '''barcode step class
    '''   
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
        self.total_num = 0
        self.clean_num = 0
        self.no_polyT_num = 0
        self.lowQual_num = 0
        self.no_linker_num = 0
        self.no_barcode_num = 0
        self.barcode_qual_Counter = Counter()
        self.umi_qual_Counter = Counter()
        self.C_U_base_Counter = Counter()
        if args.gzip:
            suffix = ".gz"
        else:
            suffix = ""
        self.out_fq2 = f'{self.outdir}/{self.sample}_2.fq{suffix}'
        self.Barcode_dict = defaultdict(int)
        self.nopolyT = args.nopolyT
        self.noLinker = args.noLinker
        self.bool_probe = False
        if args.probe_file and args.probe_file != 'None':
            self.bool_probe = True
            self.probe_count_dic = utils.genDict(dim=3)
            self.valid_count_dic = utils.genDict(dim=2)
            self.probe_dic = utils.read_fasta(args.probe_file)
            self.reads_without_probe = 0
        self.pattern = args.pattern
        self.linker = args.linker
        self.whitelist = args.whitelist
        self.lowNum = args.lowNum
        self.lowQual = args.lowQual
        self.allowNoPolyT = args.allowNoPolyT
        self.allowNoLinker = args.allowNoLinker


    def no_barcode(self, seq_arr, mis_dict, err_tolerance=1):
        tmp = [mis_dict[seq][0:2] if seq in mis_dict else (
            'X', 100) for seq in seq_arr]
        err = sum([t[1] for t in tmp])
        if err > err_tolerance:
            return True
        else:
            if err > 0:
                self.barcode_corrected_num += 1
                return ''.join([t[0] for t in tmp])
            else:
                return "correct"


    @utils.add_log
    def run(self):

        fh3 = xopen(self.out_fq2, 'w')

        if self.nopolyT:
            fh1_without_polyT = xopen(self.outdir + '/noPolyT_1.fq', 'w')
            fh2_without_polyT = xopen(self.outdir + '/noPolyT_2.fq', 'w')

        if self.noLinker:
            fh1_without_linker = xopen(self.outdir + '/noLinker_1.fq', 'w')
            fh2_without_linker = xopen(self.outdir + '/noLinker_2.fq', 'w')

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
            # generate list with mismatch 1, substitute one base in raw sequence with A,T,C,G
            if bool_whitelist:
                barcode_dict = generate_seq_dict(whitelist, n=1)
            if bool_L:
                linker_dict = generate_seq_dict(linker, n=2)
                # check linker
                check_seq(linker, pattern_dict, "L")

            fq1 = pysam.FastxFile(self.fq1_list[i], persist=False)
            fq2 = pysam.FastxFile(self.fq2_list[i], persist=False)

            for entry1 in fq1:
                entry2 = next(fq2)
                header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                if self.total_num > 0 and self.total_num % 1000000 == 0:
                    Barcode.run.logger.info(
                        f'processed reads: {utils.format_number(self.total_num)}.'
                        f'valid reads: {utils.format_number(self.clean_num)}.'
                    )

                self.total_num += 1

                # polyT filter
                if bool_T and (not self.allowNoPolyT):
                    polyT = seq_ranges(seq1, pattern_dict['T'])
                    if no_polyT(polyT):
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
                if low_qual(C_U_quals_ascii, lowQual, lowNum):
                    self.lowQual_num += 1
                    continue

                # linker filter
                barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
                raw_cb = ''.join(barcode_arr)
                if bool_L and (not self.allowNoLinker):
                    linker = seq_ranges(seq1, pattern_dict['L'])
                    if (no_linker(linker, linker_dict)):
                        self.no_linker_num += 1

                        if self.noLinker:
                            fh1_without_linker.write(
                                '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                            fh2_without_linker.write(
                                '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                        continue

                # barcode filter
                # barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
                # raw_cb = ''.join(barcode_arr)
                res = "correct"
                if bool_whitelist:
                    res = self.no_barcode(barcode_arr, barcode_dict)
                if res is True:
                    self.no_barcode_num += 1
                    continue
                elif res == "correct":
                    cb = raw_cb
                else:
                    cb = res

                umi = seq_ranges(seq1, pattern_dict['U'])
                self.Barcode_dict[cb] += 1
                self.clean_num += 1
                read_name_probe = 'None'

                if self.bool_probe:
                    # valid count
                    self.valid_count_dic[cb][umi] += 1

                    # output probe UMi and read count
                    find_probe = False
                    for probe_name in self.probe_dic:
                        probe_seq = self.probe_dic[probe_name]
                        probe_seq = probe_seq.upper()
                        if seq1.find(probe_seq) != -1:
                            self.probe_count_dic[probe_name][cb][umi] += 1
                            read_name_probe = probe_name
                            find_probe = True
                            break

                    if not find_probe:
                        self.reads_without_probe += 1

                self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])
                self.C_U_base_Counter.update(raw_cb + umi)

                # new readID: @barcode_umi_old readID
                fh3.write(f'@{cb}_{umi}_{read_name_probe}_{self.total_num}\n{seq2}\n+\n{qual2}\n')
            Barcode.run.logger.info(self.fq1_list[i] + ' finished.')
        fh3.close()

        # logging
        if self.total_num % 1000000 != 0:
            Barcode.run.logger.info(
                f'processed reads: {utils.format_number(self.total_num)}. '
                f'valid reads: {utils.format_number(self.clean_num)}. '
            )

        Barcode.run.logger.info(f'no polyT reads number : {self.no_polyT_num}')
        Barcode.run.logger.info(f'low qual reads number: {self.lowQual_num}')
        Barcode.run.logger.info(f'no_linker: {self.no_linker_num}')
        Barcode.run.logger.info(f'no_barcode: {self.no_barcode_num}')

        if self.clean_num == 0:
            raise Exception(
                'no valid reads found! please check the --chemistry parameter.')

        if self.bool_probe:
            # total probe summary
            total_umi = 0
            total_valid_read = 0
            for cb in self.valid_count_dic:
                total_umi += len(self.valid_count_dic[cb])
                total_valid_read += sum(self.valid_count_dic[cb].values())

            # probe summary
            count_list = []
            for probe_name in self.probe_dic:
                UMI_count = 0
                read_count = 0
                if probe_name in self.count_dic:
                    for cb in self.count_dic[probe_name]:
                        UMI_count += len(self.count_dic[probe_name][cb])
                        read_count += sum(self.count_dic[probe_name][cb].values())
                count_list.append(
                    {"probe_name": probe_name, "UMI_count": UMI_count, "read_count": read_count})

            df_count = pd.DataFrame(count_list, columns=[
                                    "probe_name", "read_count", "UMI_count"])

            def format_percent(x):
                x = str(round(x*100, 2))+"%"
                return x
            df_count["read_fraction"] = (
                df_count["read_count"]/total_valid_read).apply(format_percent)
            df_count["UMI_fraction"] = (
                df_count["UMI_count"]/total_umi).apply(format_percent)
            df_count.sort_values(by="UMI_count", inplace=True, ascending=False)
            df_count_file = self.outdir + '/' + self.sample + '_probe_count.tsv'
            df_count.to_csv(df_count_file, sep="\t", index=False)

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
        with open(self.outdir + '/stat.txt', 'w') as fh:
            stat_info = stat_info % (utils.format_number(self.total_num), utils.format_number(self.clean_num),
                                    cal_percent(self.clean_num), BarcodesQ30,
                                    UMIsQ30)
            stat_info = re.sub(r'^\s+', r'', stat_info, flags=re.M)
            fh.write(stat_info)
        
        self.clean_up()


    @utils.add_log
    def fastqc(self):
        cmd = ['fastqc', '-t', str(self.thread), '-o', self.outdir, self.out_fq2]
        Barcode.fastqc.logger.info('%s' % (' '.join(cmd)))
        subprocess.check_call(cmd)


@utils.add_log
def barcode(args):
    step_name = "barcode"
    barcode_obj = Barcode(args, step_name)
    barcode_obj.run()


def get_opts_barcode(parser, sub_program=True):
    parser.add_argument('--pattern', help='')
    parser.add_argument('--whitelist', help='')
    parser.add_argument('--linker', help='')
    parser.add_argument('--lowQual', type=int,
                        help='max phred of base as lowQual, default=0', default=0)
    parser.add_argument(
        '--lowNum', type=int, help='max number with lowQual allowed, default=2', default=2)
    parser.add_argument('--nopolyT', action='store_true',
                        help='output nopolyT fq')
    parser.add_argument('--noLinker', action='store_true',
                        help='output noLinker fq')
    parser.add_argument('--probe_file', help="probe fasta file")
    parser.add_argument('--allowNoPolyT', help="allow reads without polyT", action='store_true')
    parser.add_argument('--allowNoLinker', help="allow reads without correct linker", action='store_true')
    parser.add_argument('--gzip', help="output gzipped fastq", action='store_true')
    parser.add_argument(
        '--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version', default='auto')
    if sub_program:
        parser.add_argument('--fq1', help='read1 fq file', required=True)
        parser.add_argument('--fq2', help='read2 fq file', required=True)
        parser = s_common(parser)

    return parser
