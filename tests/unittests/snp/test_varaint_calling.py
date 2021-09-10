import unittest
import os
import cProfile

import celescope.snp.variant_calling as vc
from tests.unittests.__init__ import TEST_DIR_ROOT



class Test_snp(unittest.TestCase):
    def setUp(self):
        os.chdir(f'{TEST_DIR_ROOT}/snp/')
        self.sample = 'test1'
        self.vcf_file = f'{self.sample}/07.variant_calling/{self.sample}_merged.vcf'

    @unittest.skip('skip')
    def test_cell_UMI(self):
        cid = 7996
        outdir = f'{self.sample}/07.variant_calling/'
        df_UMI = vc.cell_UMI(cid, outdir, self.vcf_file)
        print(df_UMI)

    @unittest.skip('skip')
    def test_cell_UMI_2(self):
        cid = 2301
        os.chdir('/SGRNJ03/randd/RD20081701_SCOPEv2_Dynaseq/20210820/')
        sample = 'Mito_9geneMix_0812'
        outdir = f'{sample}/07.variant_calling/'
        vcf_file = f'{sample}/07.variant_calling/{sample}_merged.vcf'
        df_UMI = vc.cell_UMI(cid, outdir, vcf_file)
        print(df_UMI)


if __name__ == '__main__':
    unittest.main()