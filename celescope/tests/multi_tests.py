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
        self.test_dir = "/SGRNJ/Database/script/tests/multi_tests/"

    def _run(self, assay):
        os.chdir(self.test_dir + assay)
        print('running ' + assay)
        subprocess.check_call('sh run_shell.sh', shell=True)
        subprocess.check_call('sh sjm.sh', shell=True)
        if os.path.exists("test1"):
            shutil.rmtree("test1")
        subprocess.check_call('sh ./shell/test1.sh', shell=True)

    #@unittest.skip("pass")
    def test_rna(self):
        self._run('rna')
        assert os.path.exists("test1/06.analysis/stat.txt")

    def test_vdj(self):
        self._run('vdj')
        stat_file = glob.glob("test1/*.count_vdj/stat.txt")[0]
        assert os.path.exists(stat_file)


if __name__ == '__main__':
    unittest.main()
