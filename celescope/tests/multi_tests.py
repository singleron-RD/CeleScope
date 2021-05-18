"""
Integration tests
"""

import unittest
import os
import subprocess
import shutil
from concurrent import futures


ASSAYS = [
    'rna',
    'vdj',
    'tag',
    'capture_virus',
    #'snp',
    'fusion',
]
TEST_DIR = "/SGRNJ03/randd/user/zhouyiqi/multi_tests/"


def run_single(assay):
    """
    Returns:
        string indicates complete status
    """
    os.chdir(TEST_DIR + assay)
    print("*" * 20 + "running " + assay + "*" * 20)
    subprocess.check_call('sh run_shell.sh', shell=True)
    subprocess.check_call('sh sjm.sh', shell=True)
    if os.path.exists("test1"):
        shutil.rmtree("test1")
    try:
        subprocess.check_call('sh ./shell/test1.sh', shell=True)
    except subprocess.CalledProcessError:
        return f"{assay} failed"
    print("*" * 20 + "success " + assay + "*" * 20)
    return f"{assay} success."

class Tests(unittest.TestCase):

    def setUp(self):
        self.thread = len(ASSAYS)
        # self.thread = 1

    def test_multi(self):
        executor = futures.ProcessPoolExecutor(max_workers=self.thread)
        results = executor.map(run_single, ASSAYS)
        res_list = []
        for result in results:
            res_list.append(result)
        for result in res_list:
            print(result)
        assert not any((string.find("failed") != -1 for string in res_list))


if __name__ == '__main__':
    unittest.main()
