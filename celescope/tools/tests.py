import unittest
from collections import namedtuple

from celescope.tools.consensus import dumb_consensus, get_read_length
from celescope.tools.count import Count
from celescope.tools.step import Step


class Tests(unittest.TestCase):
    """
    Run this test under a temp folder as it will generate some files.
    """

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

    def test_get_read_length(self):
        read_list = [['AAAA', 'FFFF'], ['TTT', 'FFF'], ['CCC', 'FFF'], ['GGGGGGG', 'FFFFFFF']]
        assert get_read_length(read_list, 0.5) == 4

    def test_dumb_consensus(self):
        read_list = [('AAAA', 'FFFF'), ('TTT', 'FF;'), ('CCC', 'FFF'), ('GGGGGGG', 'FFFFFFF')]
        consensus_seq, consensus_qual, _ambiguous_base_n, _con_len = dumb_consensus(read_list, 0.5)
        print(consensus_qual)
        assert consensus_seq == 'NNNA'

    def test_correct_umi(self):
        dic = {
            "apple1": 2,
            "apple2": 30,
            "bears1": 5,
            "bears2": 10,
            "bears3": 100,
            "ccccc1": 20,
            "ccccc2": 199,
        }
        n_corrected_umi, n_corrected_read = Count.correct_umi(dic)
        sorted_dic = sorted(dic.items(), key=lambda x: x[1])
        assert sorted_dic == [('ccccc1', 20), ('apple2', 32), ('bears3', 115), ('ccccc2', 199)]
        assert n_corrected_umi == 3
        assert n_corrected_read == 2 + 5 + 10


if __name__ == '__main__':
    unittest.main()
