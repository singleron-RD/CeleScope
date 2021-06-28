"""barcode step."""

import re
from collections import Counter

import pandas as pd
import pysam
from xopen import xopen

import celescope.tools.utils as utils
from celescope.tools.__init__ import __PATTERN_DICT__
from celescope.tools.barcode import *
from celescope.tools.step import Step, s_common


class Convert(Step):

    '''
    Features

    - Demultiplex barcodes.
    - Filter invalid R1 reads, which includes:
        - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
        - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
        - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
        - Low quality reads: low sequencing quality in barcode and UMI regions.

    Output

    - `01.convert/{sample}_2.fq(.gz)`, `01.convert/{sample}_2.fq(.gz)`. Barcode and UMI are contained in the R1 reads.
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
        self.linker_corrected_num = 0
        self.total_num = 0
        self.clean_num = 0
        self.no_polyT_num = 0
        self.lowQual_num = 0
        self.no_linker_num = 0
        self.no_barcode_num = 0
        self.barcode_qual_Counter = Counter()
        self.umi_qual_Counter = Counter()
        if args.gzip:
            suffix = ".gz"
        else:
            suffix = ""
        self.out_fq1 = f'{self.outdir}/{self.sample}_1.fq{suffix}'
        self.out_fq2 = f'{self.outdir}/{self.sample}_2.fq{suffix}'
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

    @utils.add_log
    def run(self):

        outfq1 = xopen(self.out_fq1, 'w')
        outfq2 = xopen(self.out_fq2, 'w')

        if self.nopolyT:
            fh1_without_polyT = xopen(self.outdir + '/noPolyT_1.fq', 'w')
            fh2_without_polyT = xopen(self.outdir + '/noPolyT_2.fq', 'w')

        if self.noLinker:
            fh1_without_linker = xopen(self.outdir + '/noLinker_1.fq', 'w')
            fh2_without_linker = xopen(self.outdir + '/noLinker_2.fq', 'w')

        for i in range(self.fq_number):

            chemistry = self.chemistry_list[i]
            lowNum = int(self.lowNum)
            Convert.run.logger.info(f'lowQual score: {self.lowQual}')
            lowQual = int(self.lowQual)
            if chemistry == 'scopeV1':
                lowNum = min(0, lowNum)
                lowQual = max(10, lowQual)
                Convert.run.logger.info(f'scopeV1: lowNum={lowNum}, lowQual={lowQual} ')
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
                            find_probe = True
                            break

                    if not find_probe:
                        self.reads_without_probe += 1

                self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])

                outfq1.write(f'@{header1}\n{cb}{umi}\n+\n{C_U_quals_ascii}\n')

                outfq2.write(f'@{header2}\n{seq2}\n+\n{qual2}\n')

            Convert.run.logger.info(self.fq1_list[i] + ' finished.')
        outfq1.close()
        outfq2.close()

        # logging
        Convert.run.logger.info(
            f'processed reads: {utils.format_number(self.total_num)}. '
            f'valid reads: {utils.format_number(self.clean_num)}. '
        )

        Convert.run.logger.info(f'no polyT reads number : {self.no_polyT_num}')
        Convert.run.logger.info(f'low qual reads number: {self.lowQual_num}')
        Convert.run.logger.info(f'no_linker: {self.no_linker_num}')
        Convert.run.logger.info(f'no_barcode: {self.no_barcode_num}')
        Convert.run.logger.info(f'corrected linker: {self.linker_corrected_num}')
        Convert.run.logger.info(f'corrected barcode: {self.barcode_corrected_num}')

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
                if probe_name in self.probe_count_dic:
                    for cb in self.probe_count_dic[probe_name]:
                        UMI_count += len(self.probe_count_dic[probe_name][cb])
                        read_count += sum(self.probe_count_dic[probe_name][cb].values())
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
        
        # self.fastqc()
        self.clean_up()


@utils.add_log
def convert(args):
    step_name = "convert"
    convert_obj = Convert(args, step_name)
    convert_obj.run()


def get_opts_convert(parser, sub_program=True):
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

