
from celescope.tools.capture.filter import Filter, get_opts_filter


def get_opts_filter_fusion(parser, sub_program):

    get_opts_filter(parser, sub_program)

def filter_fusion(args):

    with Filter_fusion(args) as runner:
        runner.run()


class Filter_fusion(Filter):
    pass