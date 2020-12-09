from collections import defaultdict
import os
import pysam
from .utils import log, read_one_col
from .utils import seq_ranges, get_scope_bc, parse_pattern


class Chemistry():

    def __init__(self, fq1):
        self.fq1 = fq1
        self.fq1_list = fq1.split(',')
        self.nRead = 10000

    @log
    def check_chemistry(self):
        chemistry_list = []
        for fastq1 in self.fq1_list:
            print(fastq1)
            chemistry = self.get_chemistry(fastq1)
            chemistry_list.append(chemistry)
        if len(set(chemistry_list)) != 1:
            raise Exception('multiple chemistry found!' + str(chemistry_list))
        return chemistry

    @log
    def get_chemistry(self, fq1):
        '''
        'scopeV2.0.1': 'C8L16C8L16C8L1U8T18'
        'scopeV2.1.1': 'C8L16C8L16C8L1U12T18'
        'scopeV2.2.1': 'C8L16C8L16C8L1U12T18' with 4 types of linkers
        '''
        # init
        linker_4_file, _whitelist = get_scope_bc('scopeV2.2.1')
        linker_4_list, _num = read_one_col(linker_4_file)
        linker_4_dict = defaultdict(int)
        pattern_dict = parse_pattern('C8L16C8L16C8L1U12T18')
        T4_n = 0
        L57C_n = 0

        with pysam.FastxFile(fq1) as fh:
            for i in range(self.nRead):
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

        percent_T4 =  T4_n / self.nRead 
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
                raise Exception('auto chemistry detection failed!')
            elif valid_linker_type == 1:
                chemistry = 'scopeV2.1.1'
            elif valid_linker_type < 4:
                chemistry = 'scopeV2.1.1'
                Chemistry.get_chemistry.logger.warning(
                    f'chemistry scopeV2.2.1 only has {valid_linker_type} linker types!')
            else:
                chemistry = 'scopeV2.2.1'
        Chemistry.get_chemistry.logger.info(f'chemistry: {chemistry}')
        return chemistry



