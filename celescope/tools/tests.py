import unittest
import os
import pandas as pd
from collections import namedtuple
from celescope.tools.Step import Step


class Tests(unittest.TestCase):
    def setUp(self):
        pass


    def test_Step(self):
        Args = namedtuple("Args", 'sample outdir assay debug thread')
        args = Args(
            sample='test1',
            outdir='01.target_metrics',
            assay='snp',
            debug=False,
            thread=4,
        )
        step_name = "target_metrics"
        step = Step(args, step_name)

        step.add_metric(name="number value", value=1234)
        step.add_metric(name="fraction", fraction=0.123456)
        step.add_metric(name="str value", value="some str")
        step.add_metric(name="value with total", value=1234, total=3456)

        step.clean_up()


if __name__ == '__main__':
    unittest.main()