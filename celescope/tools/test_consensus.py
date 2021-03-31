import unittest
import os
from .consensus import *

class Test_Fastq(unittest.TestCase):

    #@unittest.skip('pass')
    def test_get_read_length(self):
        read_list = [['AAAA','FFFF'],['TTT','FFF'],['CCC','FFF'],['GGGGGGG','FFFFFFF']]
        assert get_read_length(read_list, 0.5) == 4
        #read_list = ['AAAA','TTT','CCC','GGGGGGG','GGGGGGG','GGGGGGG','GGGGGGG']
        #assert Fastq.get_read_length(read_list, 0.5) == 7

    #@unittest.skip('pass')
    def test_dumb_consensus(self):
        read_list = [('AAAA','FFFF'),('TTT','FF;'),('CCC','FFF'),('GGGGGGG','FFFFFFF')]
        consensus_seq, consensus_qual, ambiguous_base_n, con_len = dumb_consensus(read_list, 0.5)
        print(consensus_qual)
        assert consensus_seq == 'NNNA'

    #@unittest.skip('pass')
    def test_sorted_dumb_consensus(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/vdj/rebuild/TCR_20201221_consensus/03.consensus')
        fastq = 'TCR_20201221_consensus_sorted.fq.tmp'
        sorted_dumb_consensus(fastq, outfile = 'update_consensus.fq.gz', threshold=0.5)




if __name__ == '__main__':
    unittest.main()