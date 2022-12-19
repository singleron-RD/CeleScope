from xopen import xopen

from celescope.tools import barcode as super_barcode
from celescope.tools import utils
from celescope.tools.__init__ import PATTERN_DICT
from celescope.__init__ import HELP_DICT

import pysam 


class Barcode(super_barcode.Barcode):
    """
    ## Features
    - Demultiplex barcodes.
    - Filter invalid R2 reads, which includes:
        - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
        - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
        - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
        - Low quality reads: low sequencing quality in barcode and UMI regions.

    ## Output
    - `01.barcode/{sample}_1.fq(.gz)`. Read 1
    - `01.barcode/{sample}_2.fq(.gz)`. Dual index i5 read(Barcode only)
    - `01.barcode/{sample}_3.fq(.gz)`. Read 2
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.fq3_list = args.fq3.split(",")

    def open_files(self):

        self.out_fq3 = f'{self.out_prefix}_3.fq{self.suffix}'

        self.fh_fq1 = xopen(self.out_fq1, 'w')
        self.fh_fq2 = xopen(self.out_fq2, 'w')
        self.fh_fq3 = xopen(self.out_fq3, 'w')

        if self.nopolyT:
            self.fh_nopolyT_fq1 = xopen(self.nopolyT_1, 'w')
            self.fh_nopolyT_fq2 = xopen(self.nopolyT_2, 'w')
            self.fh_nopolyT_fq2 = xopen(self.nopolyT_3, 'w')

        if self.noLinker:
            self.fh_nolinker_fq1 = xopen(self.noLinker_1, 'w')
            self.fh_nolinker_fq2 = xopen(self.noLinker_2, 'w')
            self.fh_nolinker_fq2 = xopen(self.noLinker_3, 'w')

    def close_files(self):
        self.fh_fq1.close()
        self.fh_fq2.close()
        self.fh_fq3.close()

        if self.nopolyT:
            self.fh_nopolyT_fq1.close()
            self.fh_nopolyT_fq2.close()
            self.fh_nopolyT_fq3.close()
        
        if self.noLinker:
            self.fh_nolinker_fq1.close()
            self.fh_nolinker_fq2.close()
            self.fh_nolinker_fq3.close()

    @utils.add_log
    def add_step_metrics(self):

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

        BarcodesQ30 = sum([self.barcode_qual_Counter[k] for k in self.barcode_qual_Counter if k >= self.ord2chr(
            30)]) / float(sum(self.barcode_qual_Counter.values())) * 100
        BarcodesQ30 = round(BarcodesQ30, 2)
        BarcodesQ30_display = f'{BarcodesQ30}%'
        self.add_metric(
            name='Q30 of Barcodes',
            value=BarcodesQ30,
            display=BarcodesQ30_display,
            help_info='percent of barcode base pairs with quality scores over Q30',
        )

        try:
            UMIsQ30 = sum([self.umi_qual_Counter[k] for k in self.umi_qual_Counter if k >= self.ord2chr(
                30)]) / float(sum(self.umi_qual_Counter.values())) * 100
            UMIsQ30 = round(UMIsQ30, 2)
            UMIsQ30_display = f'{UMIsQ30}%'
    
        except ZeroDivisionError:
            UMIsQ30 = None
            UMIsQ30_display = "No UMI Pattern"

        self.add_metric(
            name='Q30 of UMIs',
            value=UMIsQ30,
            display=UMIsQ30_display,
            help_info='percent of UMI base pairs with quality scores over Q30',
        )

        self.add_metric(
            name='No PolyT Reads',
            value=self.no_polyT_num,
            total=self.total_num,
            show=False
        )

        self.add_metric(
            name='Low Quality Reads',
            value=self.lowQual_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='No Linker Reads',
            value=self.no_linker_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='No Barcode Reads',
            value=self.no_barcode_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='Corrected Linker Reads',
            value=self.linker_corrected_num,
            total=self.total_num,
            show=False,
        )

        self.add_metric(
            name='Corrected Barcode Reads',
            value=self.barcode_corrected_num,
            total=self.total_num,
            show=False,
        )

        if self.clean_num == 0:
            raise Exception('no valid reads found! please check the --chemistry parameter.' + HELP_DICT['chemistry'])
            
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

        assert "atac" in self._assay

        for i in range(self.fq_number):

            chemistry = self.chemistry_list[i]
            lowNum = int(self.lowNum)
            lowQual = int(self.lowQual)

            # get linker and whitelist
            bc_pattern = PATTERN_DICT[chemistry]
            if (bc_pattern):
                linker_file, whitelist_file = self.get_scope_bc(chemistry)
            else:
                bc_pattern = self.pattern
                linker_file = self.linker
                whitelist_file = self.whitelist
            if not bc_pattern:
                raise Exception("invalid bc_pattern!")

            pattern_dict = self.parse_pattern(bc_pattern)

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
                    pysam.FastxFile(self.fq2_list[i], persist=False) as fq2, \
                        pysam.FastxFile(self.fq3_list[i], persist=False) as fq3:


                for entry1, entry2, entry3 in zip(fq1, fq2, fq3):
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    header3, seq3, qual3 = entry3.name, entry3.sequence, entry3.quality
                    self.total_num += 1

                    # polyT filter
                    if bool_T and self.filterNoPolyT:
                        if not Barcode.check_polyT(seq2, pattern_dict):
                            self.no_polyT_num += 1
                            if self.nopolyT:
                                self.fh_nopolyT_fq1.write(
                                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                                self.fh_nopolyT_fq2.write(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                                self.fh_nopolyT_fq3.write(
                                    '@%s\n%s\n+\n%s\n' % (header3, seq3, qual3))
                            continue

                    # lowQual filter
                    C_U_quals_ascii = Barcode.get_seq_str(
                        qual2, pattern_dict['C'] + pattern_dict['U'])
                    if lowQual > 0 and Barcode.low_qual(C_U_quals_ascii, lowQual, lowNum):
                        self.lowQual_num += 1
                        continue

                    # linker filter
                    if bool_L and (not self.allowNoLinker):
                        seq_str = Barcode.get_seq_str(seq2, pattern_dict['L'])
                        bool_valid, bool_corrected, _ = Barcode.check_seq_mismatch(
                            [seq_str], linker_set_list, linker_mismatch_list)
                        if not bool_valid:
                            self.no_linker_num += 1
                            if self.noLinker:
                                self.fh_noLinker_fq1(
                                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                                self.fh_noLinker_fq2(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                                self.fh_noLinker_fq3(
                                    '@%s\n%s\n+\n%s\n' % (header3, seq3, qual3))
                            continue
                        elif bool_corrected:
                            self.linker_corrected_num += 1

                    # barcode filter
                    seq_list = self.get_seq_list(seq2, pattern_dict, 'C')

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

                    self.clean_num += 1
                    self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                    self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])

                    umi = Barcode.get_seq_str(seq2, pattern_dict['U'])
                    if umi:
                        umi += '_'

                    qual2 = len(cb) * 'F'
                    self.fh_fq3.write(f'@{cb}_{umi}{self.total_num}\n{seq3}\n+\n{qual3}\n')
                    self.fh_fq2.write(f'@{cb}_{umi}{self.total_num}\n{cb}\n+\n{qual2}\n')
                    self.fh_fq1.write(f'@{cb}_{umi}{self.total_num}\n{seq1}\n+\n{qual1}\n')                   
            
                self.run.logger.info(self.fq1_list[i] + ' finished.')

        self.close_files()
        self.add_step_metrics()


@utils.add_log
def barcode(args):
    with Barcode(args, display_title='Demultiplexing') as runner:
        runner.run()


def get_opts_barcode(parser, sub_program=True):
    super_barcode.get_opts_barcode(parser, sub_program)
    if sub_program:
        parser.add_argument('--fq3', help='R3 fastq file. Multiple files are separated by comma.', required=True)