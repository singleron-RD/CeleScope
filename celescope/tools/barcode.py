#!/bin/env python
# coding=utf8
import os
import re
import io
import gzip
import subprocess
import sys
import logging
from collections import defaultdict, Counter
from itertools import combinations, permutations, islice
from xopen import xopen
from celescope.tools.utils import format_number
from celescope.tools.report import reporter
from celescope.tools.__init__ import __PATTERN_DICT__


logger1 = logging.getLogger(__name__)
barcode_corrected_num = 0

# 定义输出格式
stat_info = '''
    Raw Reads: %s
    Valid Reads: %s(%s)
    Q30 of Barcodes: %.2f%%
    Q30 of UMIs: %.2f%%
'''


def ord2chr(q, offset=33):
    return chr(int(q) + offset)


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
                res[''.join(seq_tmp)] = (seq, mis_num, ','.join([str(i)
                                                                 for i in g]), ','.join(raw_tmp), ','.join(b))
    return(res)


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
                    logger1.warning('barcode %s, %s\n%s, %s' %
                                    (v, k, seq_dict[k], k))
                else:
                    seq_dict[k] = v
    return seq_dict


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
        logger1.error('Can not recognise pattern! %s' % pattern)
    start = 0
    for item in tmp:
        end = start + int(item[1])
        pattern_dict[item[0]].append([start, end])
        start = end
    return pattern_dict


def get_scope_bc(bctype):
    tools_path = os.path.dirname(os.path.abspath(__file__))
    root_path = os.path.dirname(tools_path)
    linker_f = os.path.join(root_path, 'data/' + bctype + '/linker_withC')
    whitelist_f = os.path.join(root_path, 'data/' + bctype + '/bclist')
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
        raise FormatError("FASTQ file ended prematurely")


def seq_ranges(seq, arr):
    # get subseq with intervals in arr and concatenate
    length = len(seq)
    for x in arr:
        if length < x[1]:
            raise Exception(f'invalid seq range {x[0]}:{x[1]} in read')
    return ''.join([seq[x[0]:x[1]]for x in arr])


def low_qual(quals, minQ='/', num=2):
    # print(ord('/')-33)           14
    return True if len([q for q in quals if q < minQ]) > num else False


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


def barcode(args):
    logger1.info('extract barcode ...!')

    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % args.outdir)

    if (args.chemistry):
        bc_pattern = __PATTERN_DICT__[args.chemistry]
        (linker, whitelist) = get_scope_bc(args.chemistry)
    else:
        bc_pattern = args.pattern
        linker = args.linker
        whitelist = args.whitelist
    if (not linker) or (not whitelist) or (not bc_pattern):
        sys.exit("invalid chemistry or [pattern,linker,whitelist]")

    # parse pattern to dict, C8L10C8L10C8U8
    # defaultdict(<type 'list'>, {'C': [[0, 8], [18, 26], [36, 44]], 'U':
    # [[44, 52]], 'L': [[8, 18], [26, 36]]})
    pattern_dict = parse_pattern(bc_pattern)
    #pattern_dict = parse_pattern(args.pattern)
    bool_T = True if 'T' in pattern_dict else False
    bool_L = True if 'L' in pattern_dict else False

    C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

    barcode_qual_Counter = Counter()
    umi_qual_Counter = Counter()
    C_U_base_Counter = Counter()
    args.lowQual = ord2chr(args.lowQual)

    # generate list with mismatch 1, substitute one base in raw sequence with
    # A,T,C,G
    barcode_dict = generate_seq_dict(whitelist, n=1)
    linker_dict = generate_seq_dict(linker, n=2)

    fq1_list = args.fq1.split(",")
    fq2_list = args.fq2.split(",")
    # merge multiple fastq files
    if len(fq1_list) > 1:
        logger1.info("merge fastq with same sample name...")
        fastq_dir = args.outdir + "/../merge_fastq"
        if not os.path.exists(fastq_dir):
            os.system('mkdir -p %s' % fastq_dir)
        fastq1_file = "{fastq_dir}/{sample}_1.fq.gz".format(
            fastq_dir=fastq_dir, sample=args.sample)
        fastq2_file = "{fastq_dir}/{sample}_2.fq.gz".format(
            fastq_dir=fastq_dir, sample=args.sample)
        fq1_cmd = "cat {fq1_files} > {fastq1_file}".format(
            fq1_files=" ".join(fq1_list), fastq1_file=fastq1_file)
        fq2_cmd = "cat {fq2_files} > {fastq2_file}".format(
            fq2_files=" ".join(fq2_list), fastq2_file=fastq2_file)
        logger1.info(fq1_cmd)
        os.system(fq1_cmd)
        logger1.info(fq2_cmd)
        os.system(fq2_cmd)
        logger1.info("merge fastq done.")
    else:
        fastq1_file = args.fq1
        fastq2_file = args.fq2

    fh1 = xopen(fastq1_file)
    fh2 = xopen(fastq2_file)
    out_fq2 = args.outdir + '/' + args.sample + '_2.fq.gz'
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

    g1 = read_fastq(fh1)
    g2 = read_fastq(fh2)
    while True:
        try:
            (header1, seq1, qual1) = next(g1)
            (header2, seq2, qual2) = next(g2)
        except BaseException:
            break

        total_num += 1
        # if total_num > 10000: total_num-= 1; break

        # polyT filter
        if bool_T:
            polyT = seq_ranges(seq1, pattern_dict['T'])
            if no_polyT(polyT):
                no_polyT_num += 1
                if args.nopolyT:
                    fh1_without_polyT.write(
                        '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                    fh2_without_polyT.write(
                        '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                continue

        # lowQual filter
        C_U_quals_ascii = seq_ranges(
            qual1, pattern_dict['C'] + pattern_dict['U'])
        # C_U_quals_ord = [ord(q) - 33 for q in C_U_quals_ascii]
        if low_qual(C_U_quals_ascii, args.lowQual, args.lowNum):
            lowQual_num += 1
            continue

        # linker filter
        barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
        raw_cb = ''.join(barcode_arr)
        if bool_L:
            linker = seq_ranges(seq1, pattern_dict['L'])
            if (no_linker(linker, linker_dict)):
                no_linker_num += 1

                if args.noLinker:
                    fh1_without_linker.write(
                        '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                    fh2_without_linker.write(
                        '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                continue

        # barcode filter
        # barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
        # raw_cb = ''.join(barcode_arr)
        res = no_barcode(barcode_arr, barcode_dict)
        if res is True:
            no_barcode_num += 1
            continue
        elif res == "correct":
            cb = raw_cb
        else:
            cb = res

        umi = seq_ranges(seq1, pattern_dict['U'])
        Barcode_dict[cb] += 1
        # new readID: @barcode_umi_old readID
        fh3.write('@{cellbarcode}_{umi}_{readID}\n{seq}\n+\n{qual}\n'.format(
            readID=header2.strip().split(' ')[0][1:], cellbarcode=cb,
            umi=umi, seq=seq2, qual=qual2))
        clean_num += 1

        barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
        umi_qual_Counter.update(C_U_quals_ascii[C_len:])
        C_U_base_Counter.update(raw_cb + umi)

    fh3.close()

    # stat
    # print(barcode_qual_Counter)
    # print(umi_qual_Counter)
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
    logger1.info('extract barcode done!')

    logger1.info('fastqc ...!')
    cmd = ['fastqc', '-t', str(args.thread), '-o', args.outdir, out_fq2]
    logger1.info('%s' % (' '.join(cmd)))
    subprocess.check_call(cmd)
    logger1.info('fastqc done!')

    t = reporter(name='barcode', assay=args.assay, sample=args.sample,
                 stat_file=args.outdir + '/stat.txt', outdir=args.outdir + '/..')
    t.get_report()


def get_opts_barcode(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument('--sample', help='sample name', required=True)
    parser.add_argument('--fq1', help='read1 fq file', required=True)
    parser.add_argument('--fq2', help='read2 fq file', required=True)
    parser.add_argument(
        '--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version')
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
    parser.add_argument('--thread', help='number of threads', default=1)
    return parser
