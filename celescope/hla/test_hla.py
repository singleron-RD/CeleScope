import os
import unittest

from celescope.hla.mapping_hla import (hla_typing, read_index, split_bam,
                                       summary)
from celescope.tools.utils import read_barcode_file


class testHLA(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/HLA/test_0908/')
        self.sample = 'pT_HLA_0812'
        self.out_bam = './/pT_HLA_0812/03.mapping_hla/bam/pT_HLA_0812.bam'
        self.match_dir = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/HLA/test_0908/match_dir/'
        self.barcodes, _ncell = read_barcode_file(self.match_dir)
        self.thread = 6
        self.mapping_outdir = f'{self.sample}/03.mapping_hla/'
        self.index_file = f'{self.mapping_outdir}/{self.sample}_cell_index.tsv'

    @unittest.skip('pass')
    def test_split_bam(self):
        outdir = self.mapping_outdir
        split_bam(self.out_bam, self.barcodes, outdir, self.sample)

    @unittest.skip('pass')
    def test_hla_typing(self):
        hla_typing(self.index_file, self.mapping_outdir, self.thread)

    @unittest.skip('pass')
    def test_read_index(self):
        read_index(self.index_file)

    # @unittest.skip('pass')
    def test_summary(self):
        summary(self.index_file, self.mapping_outdir, self.sample)


if __name__ == '__main__':
    unittest.main()
