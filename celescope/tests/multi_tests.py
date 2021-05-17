"""
Integration tests
"""

import unittest
import os
import subprocess
import shutil
import glob
from concurrent import futures


ASSAYS = {
    'rna': 'analysis',
    'vdj': 'count_vdj',
    'tag': 'count_tag',
    'capture_virus': 'analysis_capture_virus',
    'snp': 'target_metrics',
}
TEST_DIR = "/SGRNJ03/randd/user/zhouyiqi/multi_tests/"


def run_single(assay, final_stat_step):
    """
    Args:
        final_stat_step: last step with stat.txt
    """
    os.chdir(TEST_DIR + assay)
    print("*" * 20 + "running " + assay + "*" * 20)
    subprocess.check_call('sh run_shell.sh', shell=True)
    subprocess.check_call('sh sjm.sh', shell=True)
    if os.path.exists("test1"):
        shutil.rmtree("test1")
    subprocess.check_call('sh ./shell/test1.sh', shell=True)
    stat_file = glob.glob(f"test1/*.{final_stat_step}/stat.txt")[0]
    assert os.path.exists(stat_file)
    print("*" * 20 + "done " + assay + "*" * 20)
    return assay

class Tests(unittest.TestCase):

    def setUp(self):
        pass

    def test_multi(self):
        thread = len(ASSAYS)
        executor = futures.ProcessPoolExecutor(max_workers=thread)
        results = executor.map(run_single, ASSAYS.keys(), ASSAYS.values())
        for result in results:
            print(result + ' done')


if __name__ == '__main__':
    unittest.main()
