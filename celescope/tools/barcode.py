import os
import re
import io
import gzip
import subprocess
import sys
import glob
import pandas as pd
from collections import defaultdict, Counter
from itertools import combinations, permutations, islice
from xopen import xopen
from celescope.tools.utils import *
from celescope.tools.report import reporter
from celescope.tools.__init__ import __PATTERN_DICT__
from .Chemistry import Chemistry


barcode_corrected_num = 0

# 定义输出格式
stat_info = '''
    Raw Reads: %s
    Valid Reads: %s(%s)
    Q30 of Barcodes: %.2f%%
    Q30 of UMIs: %.2f%%
'''


# 生成错配字典
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


@add_log
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


@add_log
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


def no_barcode(seq_arr, mis_dict, err_tolerance=1):
    global barcode_corrected_num
    tmp = [mis_dict[seq][0:2] if seq in mis_dict else (
        'X', 100) for seq in seq_arr]
    err = sum([t[1] for t in tmp])
    if err > err_tolerance:
        return True
    else:
        if err > 0:
            barcode_corrected_num += 1
            return ''.join([t[0] for t in tmp])
        else:
            return "correct"


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


@add_log
def merge_fastq(fq1, fq2, sample, outdir):
    '''
    merge fastq if len(fq1) > 1
    '''
    fq1_list = fq1.split(",")
    fq2_list = fq2.split(",")
    if len(fq1_list) == 0:
        raise Exception('empty fastq file path!')
    elif len(fq1_list) == 1:
        fq1_file, fq2_file = fq1, fq2
    else:
        fastq_dir = f'{outdir}/../merge_fastq'
        if not os.path.exists(fastq_dir):
            os.system('mkdir -p %s' % fastq_dir)
        fq1_file = f"{fastq_dir}/{sample}_1.fq.gz"
        fq2_file = f"{fastq_dir}/{sample}_2.fq.gz"
        fq1_files = " ".join(fq1_list)
        fq2_files = " ".join(fq2_list)
        fq1_cmd = f"cat {fq1_files} > {fq1_file}"
        fq2_cmd = f"cat {fq2_files} > {fq2_file}"
        merge_fastq.logger.info(fq1_cmd)
        os.system(fq1_cmd)
        merge_fastq.logger.info(fq2_cmd)
        os.system(fq2_cmd)
    return fq1_file, fq2_file


def polyT_filter(pattern_dict, args, header1, seq1, qual1, header2, seq2, qual2, bool_T, bool_L,
                 lowQual, lowNum, linker_dict, bool_whitelist, barcode_dict):

    if bool_T and (not args.allowNoPolyT):
        polyT = seq_ranges(seq1, pattern_dict['T'])

        if no_polyT(polyT):
            if args.nopolyT:
                fh1_without_polyT = xopen(args.outdir + '/noPolyT_1.fq', 'a')
                fh2_without_polyT = xopen(args.outdir + '/noPolyT_2.fq', 'a')
                fh1_without_polyT.write(
                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                fh2_without_polyT.write(
                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
            return True
        else:
            return False
    else:
        return False


def lowQual_filter(pattern_dict, args, header1, seq1, qual1, header2, seq2, qual2, bool_T, bool_L,
                   lowQual, lowNum, linker_dict, bool_whitelist, barcode_dict):
    # lowQual filter
    C_U_quals_ascii = seq_ranges(
        qual1, pattern_dict['C'] + pattern_dict['U'])
    # C_U_quals_ord = [ord(q) - 33 for q in C_U_quals_ascii]
    if low_qual(C_U_quals_ascii, lowQual, lowNum):
        return True, C_U_quals_ascii, "Qual"
    else:
        return False, C_U_quals_ascii, "Qual"


def linker_filter(pattern_dict, args, header1, seq1, qual1, header2, seq2, qual2, bool_T, bool_L,
                  lowQual, lowNum, linker_dict, bool_whitelist, barcode_dict):

    # linker filter
    if bool_L and (not args.allowNoLinker):
        linker = seq_ranges(seq1, pattern_dict['L'])

        if (no_linker(linker, linker_dict)):
            if args.noLinker:
                fh1_without_linker = xopen(args.outdir + '/noLinker_1.fq', 'a')
                fh2_without_linker = xopen(args.outdir + '/noLinker_2.fq', 'a')
                fh1_without_linker.write(
                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                fh2_without_linker.write(
                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
            return True, "linker"
        else:
            return False, "linker"
    else:
        return False, 'linker'


def barcode_filter(pattern_dict, args, header1, seq1, qual1, header2, seq2, qual2, bool_T, bool_L,
                   lowQual, lowNum, linker_dict, bool_whitelist, barcode_dict):
    # barcode filter
    # barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
    # raw_cb = ''.join(barcode_arr)
    barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
    raw_cb = ''.join(barcode_arr)
    res = "correct"
    if bool_whitelist:
        res = no_barcode(barcode_arr, barcode_dict)
    if res is True:
        return True, "barcode"
    elif res == "correct":
        cb = raw_cb
        return False, cb, raw_cb, "barcode"
    else:
        cb = res
        return False, cb, raw_cb, "barcode"


@add_log
def barcode(args):
    # init
    outdir = args.outdir
    sample = args.sample
    fq1 = args.fq1
    fq2 = args.fq2
    fq1_list = fq1.split(",")
    fq2_list = fq2.split(",")
    filter_order = args.filter_order
    filter_order = filter_order.split(",")
    fq_number = len(fq1_list)

    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % args.outdir)

    # get chemistry
    if args.chemistry == 'auto':
        ch = Chemistry(fq1)
        chemistry_list = ch.check_chemistry()
    else:
        chemistry_list = [args.chemistry] * fq_number

    barcode_qual_Counter = Counter()
    umi_qual_Counter = Counter()
    C_U_base_Counter = Counter()

    # prepare
    if not args.not_gzip:
        suffix = ".gz"
    else:
        suffix = ""
    out_fq2 = f'{args.outdir}/{args.sample}_2.fq{suffix}'
    fh3 = xopen(out_fq2, 'w')

    (total_num, clean_num, no_polyT_num, lowQual_num,
        no_linker_num, no_barcode_num) = (0, 0, 0, 0, 0, 0)
    Barcode_dict = defaultdict(int)

    if args.nopolyT:
        fh1_without_polyT = xopen(args.outdir + '/noPolyT_1.fq', 'w')
        fh2_without_polyT = xopen(args.outdir + '/noPolyT_2.fq', 'w')

    if args.noLinker:
        fh1_without_linker = xopen(args.outdir + '/noLinker_1.fq', 'w')
        fh2_without_linker = xopen(args.outdir + '/noLinker_2.fq', 'w')

    bool_probe = False
    if args.probe_file and args.probe_file != 'None':
        bool_probe = True
        count_dic = genDict(dim=3)
        valid_count_dic = genDict(dim=2)
        probe_dic = read_fasta(args.probe_file)
        reads_without_probe = 0

    filter_dict = {'polyT': [polyT_filter, no_polyT_num],
                   'Qual': [lowQual_filter, lowQual_num],
                   'linker': [linker_filter, no_linker_num],
                   'barcode': [barcode_filter, no_barcode_num]}

    # process
    if fq_number != len(fq2_list):
        raise Exception('fastq1 and fastq2 do not have same file number!')
    for i in range(fq_number):

        chemistry = chemistry_list[i]
        lowNum = int(args.lowNum)
        print(args.lowQual)
        lowQual = int(args.lowQual)
        if chemistry == 'scopeV1':
            lowNum = min(0, lowNum)
            lowQual = max(10, lowQual)
            barcode.logger.info(f'scopeV1: lowNum={lowNum}, lowQual={lowQual} ')
        # get linker and whitelist
        bc_pattern = __PATTERN_DICT__[chemistry]
        if bc_pattern:
            (linker, whitelist) = get_scope_bc(chemistry)
        else:
            bc_pattern = args.pattern
            linker = args.linker
            whitelist = args.whitelist
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

        fh1 = xopen(fq1_list[i])
        fh2 = xopen(fq2_list[i])
        g1 = read_fastq(fh1)
        g2 = read_fastq(fh2)

        while True:
            try:
                (header1, seq1, qual1) = next(g1)
                (header2, seq2, qual2) = next(g2)
            except BaseException:
                break
            if total_num > 0 and total_num % 1000000 == 0:
                barcode.logger.info(
                    f'processed reads: {format_number(total_num)}.'
                    f'valid reads: {format_number(clean_num)}.'
                )

            total_num += 1

            # polyT filter
            filter1 = filter_dict[filter_order[0]][0](pattern_dict, args, header1, seq1, qual1, header2, seq2, qual2,
                                                      bool_T, bool_L, lowQual, lowNum, linker_dict,
                                                      bool_whitelist, barcode_dict)

            if filter_order[0] == 'polyT':
                if filter1:
                    filter_dict[filter_order[0]][1] += 1
                    continue
                else:
                    pass
            else:
                if filter1[0]:
                    filter_dict[filter_order[0]][1] += 1
                    continue
                else:
                    pass

            if filter_order[0] == 'barcode':
                cb = filter1[1]
                raw_cb = filter1[2]
            if filter_order[0] == 'Qual':
                C_U_quals_ascii = filter1[1]

            # lowQual filter
            filter2 = filter_dict[filter_order[1]][0](pattern_dict, args, header1, seq1, qual1, header2, seq2, qual2,
                                                      bool_T, bool_L, lowQual, lowNum, linker_dict,
                                                      bool_whitelist, barcode_dict)
            if filter_order[1] == 'polyT':
                if filter2:
                    filter_dict[filter_order[1]][1] += 1
                    continue
                else:
                    pass
            else:
                if filter2[0]:
                    filter_dict[filter_order[1]][1] += 1
                    continue
                else:
                    pass

            if filter_order[1] == 'barcode':
                cb = filter2[1]
                raw_cb = filter2[2]
            if filter_order[1] == 'Qual':
                C_U_quals_ascii = filter2[1]

            # linker filter
            filter3 = filter_dict[filter_order[2]][0](pattern_dict, args, header1, seq1, qual1, header2, seq2, qual2,
                                                      bool_T, bool_L, lowQual, lowNum, linker_dict,
                                                      bool_whitelist, barcode_dict)
            if filter_order[2] == 'polyT':
                if filter3:
                    filter_dict[filter_order[2]][1] += 1
                    continue
                else:
                    pass
            else:
                if filter3[0]:
                    filter_dict[filter_order[2]][1] += 1
                    continue
                else:
                    pass

            if filter_order[2] == 'barcode':
                cb = filter3[1]
                raw_cb = filter3[2]
            if filter_order[2] == 'Qual':
                C_U_quals_ascii = filter3[1]

            # barcode filter
            # barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
            # raw_cb = ''.join(barcode_arr)
            filter4 = filter_dict[filter_order[3]][0](pattern_dict, args, header1, seq1, qual1, header2, seq2, qual2,
                                                      bool_T, bool_L, lowQual, lowNum, linker_dict,
                                                      bool_whitelist, barcode_dict)
            if filter_order[3] == 'polyT':
                if filter4:
                    filter_dict[filter_order[3]][1] += 1
                    continue
                else:
                    pass
            else:
                if filter4[0]:
                    filter_dict[filter_order[3]][1] += 1
                    continue
                else:
                    pass

            if filter_order[3] == 'barcode':
                cb = filter4[1]
                raw_cb = filter4[2]
            if filter_order[3] == 'Qual':
                C_U_quals_ascii = filter4[1]

            umi = seq_ranges(seq1, pattern_dict['U'])
            Barcode_dict[cb] += 1
            clean_num += 1
            read_name_probe = 'None'

            if bool_probe:
                # valid count
                valid_count_dic[cb][umi] += 1

                # output probe UMi and read count
                find_probe = False
                for probe_name in probe_dic:
                    probe_seq = probe_dic[probe_name]
                    probe_seq = probe_seq.upper()
                    if seq1.find(probe_seq) != -1:
                        count_dic[probe_name][cb][umi] += 1
                        read_name_probe = probe_name
                        find_probe = True
                        break

                if not find_probe:
                    reads_without_probe += 1
            barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
            umi_qual_Counter.update(C_U_quals_ascii[C_len:])
            C_U_base_Counter.update(raw_cb + umi)

            # new readID: @barcode_umi_old readID
            fh3.write(f'@{cb}_{umi}_{read_name_probe}_{total_num}\n{seq2}\n+\n{qual2}\n')
        barcode.logger.info(fq1_list[i] + ' finished.')

    fh3.close()

    # logging
    if total_num % 1000000 != 0:
        barcode.logger.info(
            f'processed reads: {format_number(total_num)}. '
            f'valid reads: {format_number(clean_num)}. '
        )
    filter_count_dict = {}
    for i in filter_order:
        filter_count_dict[i] = filter_dict[i][1]
    no_linker_num = filter_count_dict['linker']
    no_polyT_num = filter_count_dict['polyT']
    lowQual_num = filter_count_dict['Qual']
    no_barcode_num = filter_count_dict['barcode']
    barcode.logger.info(f'no polyT reads number : {no_polyT_num}')
    barcode.logger.info(f'low qual reads number: {lowQual_num}')
    barcode.logger.info(f'no_linker: {no_linker_num}')
    barcode.logger.info(f'no_barcode: {no_barcode_num}')

    if clean_num == 0:
        raise Exception(
            'no valid reads found! please check the --chemistry parameter.')

    if bool_probe:
        # total probe summary
        total_umi = 0
        total_valid_read = 0
        for cb in valid_count_dic:
            total_umi += len(valid_count_dic[cb])
            total_valid_read += sum(valid_count_dic[cb].values())
        barcode.logger.info("total umi:"+str(total_umi))
        barcode.logger.info("total valid read:"+str(total_valid_read))
        barcode.logger.info("reads without probe:"+str(reads_without_probe))

        # probe summary
        count_list = []
        for probe_name in probe_dic:
            UMI_count = 0
            read_count = 0
            if probe_name in count_dic:
                for cb in count_dic[probe_name]:
                    UMI_count += len(count_dic[probe_name][cb])
                    read_count += sum(count_dic[probe_name][cb].values())
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
        df_count_file = args.outdir + '/' + args.sample + '_probe_count.tsv'
        df_count.to_csv(df_count_file, sep="\t", index=False)

    # stat
    BarcodesQ30 = sum([barcode_qual_Counter[k] for k in barcode_qual_Counter if k >= ord2chr(
        30)]) / float(sum(barcode_qual_Counter.values())) * 100
    UMIsQ30 = sum([umi_qual_Counter[k] for k in umi_qual_Counter if k >= ord2chr(
        30)]) / float(sum(umi_qual_Counter.values())) * 100

    global stat_info
    def cal_percent(x): return "{:.2%}".format((x + 0.0) / total_num)
    with open(args.outdir + '/stat.txt', 'w') as fh:
        """
        Raw Reads: %s
        Valid Reads: %s(%s)
        Q30 of Barcodes: %.2f%%
        Q30 of UMIs: %.2f%%
        """
        stat_info = stat_info % (format_number(total_num), format_number(clean_num),
                                 cal_percent(clean_num), BarcodesQ30,
                                 UMIsQ30)
        stat_info = re.sub(r'^\s+', r'', stat_info, flags=re.M)
        fh.write(stat_info)

    barcode.logger.info('fastqc ...!')
    cmd = ['fastqc', '-t', str(args.thread), '-o', args.outdir, out_fq2]
    barcode.logger.info('%s' % (' '.join(cmd)))
    subprocess.check_call(cmd)
    barcode.logger.info('fastqc done!')

    t = reporter(name='barcode', assay=args.assay, sample=args.sample,
                 stat_file=args.outdir + '/stat.txt', outdir=args.outdir + '/..')
    t.get_report()


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
    parser.add_argument('--not_gzip', help="output fastq without gzip", action='store_true')
    parser.add_argument('--filter_order', help='filter order', default='polyT,Qual,linker,barcode')
    if sub_program:
        parser.add_argument('--fq1', help='read1 fq file', required=True)
        parser.add_argument('--fq2', help='read2 fq file', required=True)
        parser.add_argument(
            '--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version', default='auto')
        parser = s_common(parser)

    return parser
