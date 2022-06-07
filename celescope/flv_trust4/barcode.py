from celescope.tools.barcode import *
from celescope.tools import barcode as tools_barcode


class Barcode(Step):
    """
    ## Features

    - Demultiplex barcodes.
    - Filter invalid R1 reads, which includes:
        - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
        - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
        - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
        - Low quality reads: low sequencing quality in barcode and UMI regions.
        - Write coverage is capped to a maximum of 80,000 reads per barcode.
        - Include only reads from match barcodes.
        - Reverse complement the barcode to match with sc-RNA.

    ## Output

    - `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
    the read name is `{barcode}_{UMI}_{read ID}`.
    - `01.barcode/{sample}_1.fq(.gz)` Write barcode and umi to R1 read(can be directly used as input file of TRUST4).
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

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
        # check flv
        for chemistry in self.chemistry_list:
            if chemistry != 'flv':
                raise Exception('chemistry should be `flv`!')

        self.barcode_corrected_num = 0
        self.linker_corrected_num = 0
        self.total_num = 0
        self.clean_num = 0
        self.match_num = 0
        self.no_polyT_num = 0
        self.lowQual_num = 0
        self.no_linker_num = 0
        self.no_barcode_num = 0
        self.barcode_qual_Counter = Counter()
        self.umi_qual_Counter = Counter()
        self.barcode_read_Counter = Counter()
        self.pattern = args.pattern
        self.linker = args.linker
        self.whitelist = args.whitelist
        self.lowNum = args.lowNum
        self.lowQual = args.lowQual
        self.allowNoPolyT = args.allowNoPolyT
        self.allowNoLinker = args.allowNoLinker
        self.nopolyT = args.nopolyT  # true == output nopolyT reads
        self.noLinker = args.noLinker
        self.match_barcodes, _ = utils.get_barcode_from_match_dir(args.match_dir)

        # out file
        if args.gzip:
            suffix = ".gz"
        else:
            suffix = ""
        self.out_fq2 = f'{self.out_prefix}_2.fq{suffix}'
        self.out_fq1 = f'{self.out_prefix}_1.fq{suffix}'
        if self.nopolyT:
            self.nopolyT_1 = f'{self.out_prefix}_noPolyT_1.fq'
            self.nopolyT_2 = f'{self.out_prefix}_noPolyT_2.fq'
        if self.noLinker:
            self.noLinker_1 = f'{self.out_prefix}_noLinker_1.fq'
            self.noLinker_2 = f'{self.out_prefix}_noLinker_2.fq'

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

        out_fq1 = xopen(self.out_fq1, 'w')
        out_fq2 = xopen(self.out_fq2, 'w')


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
            bc_pattern = PATTERN_DICT[chemistry]
            if (bc_pattern):
                linker_file, whitelist_file = get_scope_bc(chemistry)
            else:
                bc_pattern = self.pattern
                linker_file = self.linker
                whitelist_file = self.whitelist
            if not bc_pattern:
                raise Exception("invalid bc_pattern!")

            pattern_dict = parse_pattern(bc_pattern)

            bool_T = True if 'T' in pattern_dict else False
            bool_L = True if 'L' in pattern_dict else False
            bool_whitelist = (whitelist_file is not None) and whitelist_file != "None"
            C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

            if bool_whitelist:
                barcode_set_list, barcode_mismatch_list = parse_whitelist_file(whitelist_file,
                                                                               n_mismatch=1, n_repeat=len(pattern_dict['C']))
            if bool_L:
                linker_set_list, linker_mismatch_list = parse_linker_file(linker_file)

            with pysam.FastxFile(self.fq1_list[i], persist=False) as fq1, \
                    pysam.FastxFile(self.fq2_list[i], persist=False) as fq2:
                for entry1, entry2 in zip(fq1, fq2):
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    self.total_num += 1

                    # polyT filter
                    if bool_T and (not self.allowNoPolyT):
                        polyT = get_seq_str(seq1, pattern_dict['T'])
                        if polyT.count('T') < MIN_T:
                            self.no_polyT_num += 1
                            if self.nopolyT:
                                fh1_without_polyT.write(
                                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                                fh2_without_polyT.write(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue

                    # lowQual filter
                    C_U_quals_ascii = get_seq_str(
                        qual1, pattern_dict['C'] + pattern_dict['U'])
                    # C_U_quals_ord = [ord(q) - 33 for q in C_U_quals_ascii]
                    if lowQual > 0 and low_qual(C_U_quals_ascii, lowQual, lowNum):
                        self.lowQual_num += 1
                        continue

                    # linker filter
                    if bool_L and (not self.allowNoLinker):
                        seq_str = get_seq_str(seq1, pattern_dict['L'])
                        bool_valid, bool_corrected, _ = check_seq_mismatch(
                            [seq_str], linker_set_list, linker_mismatch_list)
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
                    # flv barcode is reverse complement of the flv_rna barcode
                    if bool_whitelist:
                        bool_valid, bool_corrected, corrected_seq = check_seq_mismatch(
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
                        if self.barcode_read_Counter[cb] <= 80000:
                            umi = get_seq_str(seq1, pattern_dict['U'])
                            qual = 'F' * len(cb + umi)
                            out_fq2.write(f'@{cb}_{umi}_{self.total_num}\n{seq2}\n+\n{qual2}\n')
                            out_fq1.write(f'@{cb}_{umi}_{self.total_num}\n{cb}{umi}\n+\n{qual}\n')
                        
            self.run.logger.info(self.fq1_list[i] + ' finished.')
        out_fq2.close()
        out_fq1.close()

        # logging
        self.run.logger.info(
            f'processed reads: {utils.format_number(self.total_num)}. '
            f'valid reads: {utils.format_number(self.clean_num)}. '
        )

        self.run.logger.info(f'no polyT reads number : {self.no_polyT_num}')
        self.run.logger.info(f'low qual reads number: {self.lowQual_num}')
        self.run.logger.info(f'no_linker: {self.no_linker_num}')
        self.run.logger.info(f'no_barcode: {self.no_barcode_num}')
        self.run.logger.info(f'corrected linker: {self.linker_corrected_num}')
        self.run.logger.info(f'corrected barcode: {self.barcode_corrected_num}')

        if self.clean_num == 0:
            raise Exception(
                'no valid reads found! please check the --chemistry parameter.')

        # stat
        BarcodesQ30 = sum([self.barcode_qual_Counter[k] for k in self.barcode_qual_Counter if k >= ord2chr(
            30)]) / float(sum(self.barcode_qual_Counter.values())) * 100
        BarcodesQ30 = round(BarcodesQ30, 2)
        BarcodesQ30_display = f'{BarcodesQ30}%'
        UMIsQ30 = sum([self.umi_qual_Counter[k] for k in self.umi_qual_Counter if k >= ord2chr(
            30)]) / float(sum(self.umi_qual_Counter.values())) * 100
        UMIsQ30 = round(UMIsQ30, 2)
        UMIsQ30_display = f'{UMIsQ30}%'

        self.add_metric(
            name='Raw Reads',
            value=self.total_num,
            help_info='total reads from FASTQ files'
        )
        self.add_metric(
            name='Valid Reads',
            value=self.clean_num,
            total=self.total_num,
            help_info='reads pass filtering(filtered: reads without poly T, reads without linker, reads without correct barcode or low quality reads)'
        )
        self.add_metric(
            name='Valid Matched Reads',
            value=self.match_num,
            total=self.total_num,
            help_info='reads match with flv_rna cell barcodes'
        )
        self.add_metric(
            name='Q30 of Barcodes',
            value=BarcodesQ30,
            display=BarcodesQ30_display,
            help_info='percent of barcode base pairs with quality scores over Q30',
        )
        self.add_metric(
            name='Q30 of UMIs',
            value=UMIsQ30,
            display=UMIsQ30_display,
            help_info='percent of UMI base pairs with quality scores over Q30',
        )


@utils.add_log
def barcode(args):
    with Barcode(args, display_title='Demultiplexing') as runner:
        runner.run()


def get_opts_barcode(parser, sub_program=True):
    tools_barcode.get_opts_barcode(parser, sub_program)
    if sub_program:
        parser.add_argument('--match_dir', help='Matched scRNA-seq directory.', required=True)
