import unittest
import os
import pandas as pd
from celescope.hla.mapping_hla import split_bam, hla_typing, summary
from celescope.tools.utils import read_barcode_file


class testHLA(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/HLA/sjm_0903')
        self.sample = 'pT_HLA_0812'
        self.out_bam = './/pT_HLA_0812/03.mapping_hla/bam/pT_HLA_0812.bam'
        #self.match_dir = '/SGRNJ02/RandD4/RD2019016/20200804/PBMC_pTLib/'
        self.match_dir = '/SGRNJ02/RandD4/RD2019016/20200811/pT_0730Lib'
        self.barcodes, _ncell = read_barcode_file(self.match_dir)
        self.thread = 6
        self.mapping_outdir = f'{self.sample}/03.mapping_hla/'

    @unittest.skip('pass')
    def test_split_bam(self):
        outdir = self.mapping_outdir
        split_bam(self.out_bam, outdir, self.barcodes, self.sample)

    @unittest.skip('pass')
    def test_hla_typing(self):
        outdir = self.mapping_outdir
        index_file = f'{outdir}/{self.sample}_cell_index.tsv'
        index_df = pd.read_csv(index_file, sep='\t', index_col=0)
        index_df = index_df.T
        index_dict = index_df.to_dict()
        hla_typing(index_dict, outdir, self.thread)

    def test_summary(self):
        outdir = self.mapping_outdir
        index_file = f'{outdir}/{self.sample}_cell_index.tsv'
        index_df = pd.read_csv(index_file, sep='\t', index_col=0)
        index_df = index_df.T
        index_dict = index_df.to_dict()
        valid_index_dict = {}
        for index in index_dict:
            if index_dict[index]['valid']:
                valid_index_dict[index] = index_dict[index]['barcode']
        summary(valid_index_dict, outdir, self.sample)
        


if __name__ == '__main__':
    unittest.main()
