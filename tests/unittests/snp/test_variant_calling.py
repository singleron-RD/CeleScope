import unittest
import os
from collections import namedtuple
from celescope.snp.variant_calling import Variant_calling

ROOT_DIR = os.path.dirname(__file__)


class Test_variant_calling(unittest.TestCase):
    def setUp(self):
        os.chdir(ROOT_DIR)
        Args = namedtuple("Args", "thread outdir sample assay debug " + "genomeDir vcf bam match_dir")
        self.args = Args(
            thread=10,
            outdir="./test_output/07.variant_calling",
            sample="test1",
            assay="snp",
            debug=False,
            genomeDir="/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92",
            vcf=None,
            bam="./test_data/06.target_metrics/subset_filter.bam",
            match_dir="./test_data/match_dir",
        )

    def test_run(self):
        obj = Variant_calling(self.args, "variant_calling")
        '''
        obj.SplitNCigarReads()
        obj.split_bam()
        obj.call_all_snp()
        if obj.vcf_bool:
            obj.add_VID()
        else:
            obj.merge_vcf()
        '''
        obj.write_VID_file()
        obj.get_UMI()
        obj.write_support_matrix()
        obj._clean_up()

