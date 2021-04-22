"""
Integration tests
"""
import unittest
import os
import subprocess
import shutil
import glob



class Tests(unittest.TestCase):

    def setUp(self):
        self.test_dir = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/multi_tests/"

    #@unittest.skip("pass")
    def test_rna(self):
        assay = 'rna'
        os.chdir(self.test_dir + assay)
        subprocess.check_call('sh run_shell.sh', shell=True)
        if os.path.exists("test1"):
            shutil.rmtree("test1")
        subprocess.check_call('sh ./shell/test1.sh', shell=True)
        assert os.path.exists("test1/06.analysis/stat.txt")

    def test_vdj(self):
        assay = 'vdj'
        os.chdir(self.test_dir + assay)
        subprocess.check_call('sh run_shell.sh', shell=True)
        if os.path.exists("test1"):
            shutil.rmtree("test1")
        subprocess.check_call('sh ./shell/test1.sh', shell=True)
        stat_file = glob.glob("test1/*.count_vdj/stat.txt")[0]


if __name__ == '__main__':
    unittest.main()
