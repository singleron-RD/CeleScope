
from celescope.tools.Fastq import Fastq
import unittest
import os

class Test_Fastq(unittest.TestCase):

    #@unittest.skip('pass')
    def test_get_read_length(self):
        read_list = ['AAAA','TTT','CCC','GGGGGGG']
        assert Fastq.get_read_length(read_list, 0.5) == 4
        read_list = ['AAAA','TTT','CCC','GGGGGGG','GGGGGGG','GGGGGGG','GGGGGGG']
        assert Fastq.get_read_length(read_list, 0.5) == 7

    def test_dumb_consensus(self):
        read_list = ['AAAA','TTT','CCC','GGGGGGG']
        print (Fastq.dumb_consensus(read_list, 0.5))
        assert Fastq.dumb_consensus(read_list, 0.5) == 'NNNA'
        read_list = ['AAAA','TTT','CCC','GGGGGGG','GGGGGGG','GGGGGGG','GGGGGAG','ATCGATCG']
        print (Fastq.dumb_consensus(read_list, 0.5))
        assert Fastq.dumb_consensus(read_list, 0.5) == 'GGGGGGG'

    def test_write_consensus_fasta(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/vdj/TCR/')
        outdir = 'TCR_sub/03.mapping_vdj'
        sample = 'TCR_sub'
        fq = 'TCR_sub/02.cutadapt/TCR_sub_clean_2.fq.gz'
        fq_obj = Fastq(fq)
        fq_obj.umi_dumb_consensus()
        fq_obj.write_consensus_fasta(outdir,sample)



if __name__ == '__main__':
    unittest.main()