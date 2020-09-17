import unittest
import os
import pandas as pd
from celescope.fusion.count_fusion import read_pos


class testHLA(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/CeleScope-Fusion/sjm_0818')
        self.fusion_pos = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/CeleScope-Fusion/pipe_0722/fusion_pos.txt'

    def test_read_pos(self):
        print(read_pos(self.fusion_pos))


if __name__ == '__main__':
    unittest.main()
