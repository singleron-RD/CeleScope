
from celescope.tools.Fastq import Fastq
import unittest
import os

class Test_Fastq(unittest.TestCase):

    #@unittest.skip('pass')
    def test_get_read_length(self):
        read_list = [['AAAA','FFFF'],['TTT','FFF'],['CCC','FFF'],['GGGGGGG','FFFFFFF']]
        assert Fastq.get_read_length(read_list, 0.5) == 4
        #read_list = ['AAAA','TTT','CCC','GGGGGGG','GGGGGGG','GGGGGGG','GGGGGGG']
        #assert Fastq.get_read_length(read_list, 0.5) == 7

    @unittest.skip('pass')
    def test_dumb_consensus(self):
        read_list = ['AAAA','TTT','CCC','GGGGGGG']
        print (Fastq.dumb_consensus(read_list, 0.5))
        assert Fastq.dumb_consensus(read_list, 0.5) == 'NNNA'
        read_list = ['AAAA','TTT','CCC','GGGGGGG','GGGGGGG','GGGGGGG','GGGGGAG','ATCGATCG']
        print (Fastq.dumb_consensus(read_list, 0.5))
        assert Fastq.dumb_consensus(read_list, 0.5) == 'GGGGGGG'

    @unittest.skip('pass')
    def test_write_consensus_fasta(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/vdj/TCR/')
        outdir = 'TCR_sub/03.mapping_vdj'
        sample = 'TCR_sub'
        fq = 'TCR_sub/02.cutadapt/TCR_sub_clean_2.fq.gz'
        fq_obj = Fastq(fq)
        fq_obj.umi_dumb_consensus()
        fq_obj.write_consensus_fasta(outdir,sample)
        fq_obj.write_consensus_fastq(outdir,sample)

    @unittest.skip('pass')
    def test_concurrent(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/variant/20201224_indel/consensus/')
        outdir = './S32_SUR0528_drug_S203_TS/02.cutadapt'
        sample = 'S32_SUR0528_drug_S203_TS'
        fq = './S32_SUR0528_drug_S203_TS/02.cutadapt/test_1M.fq'
        fq_obj = Fastq(fq)
        fq_obj.umi_dumb_consensus_concurrent(thread=5)
        #fq_obj.write_consensus_fasta(outdir,sample)
        #fq_obj.write_consensus_fastq(outdir,sample)

    @unittest.skip('pass')
    def test_single(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/variant/20201224_indel/consensus/')
        outdir = './S32_SUR0528_drug_S203_TS/02.cutadapt'
        sample = 'S32_SUR0528_drug_S203_TS'
        fq = './S32_SUR0528_drug_S203_TS/02.cutadapt/test_1M.fq'
        fq_obj = Fastq(fq)
        fq_obj.umi_dumb_consensus()
        #fq_obj.write_consensus_fasta(outdir,sample)
        #fq_obj.write_consensus_fastq(outdir,sample)

    def test_sorted_dumb_consensus(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/test1/01.barcode')
        fastq = 'sorted.fastq'
        fq_css = Fastq(fastq)
        fq_css.sorted_dumb_consensus(outfile = 'consensus.fq.gz')




if __name__ == '__main__':
    unittest.main()