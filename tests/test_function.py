import os
import unittest
from collections import namedtuple

from celescope.tools.step import Step


class Tests(unittest.TestCase):

    def setUp(self):
        pass

    @unittest.skip("tested")
    def test_stat_to_metric(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/multi_tests/rna')
        args_dict = {
            'sample': 'test1',
            'assay': 'rna',
            'thread': 1,
            'outdir': 'test1/06.analysis',
            'debug': True,
        }
        Args = namedtuple('Args', list(args_dict.keys()))
        args = Args(**args_dict)

        obj = Step(args, 'analysis')
        obj.stat_to_metric()
        print(obj.__content_dict['metric'])

    def test_test(self):
        assert 0 == 0
