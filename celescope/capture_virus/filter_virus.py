
from celescope.tools.capture.filter import Filter, get_opts_filter


def get_opts_filter_virus(parser, sub_program):

    get_opts_filter(parser, sub_program)

def filter_virus(args):

    with Filter_virus(args, display_title='Filtering') as runner:
        runner.run()


class Filter_virus(Filter):
    pass
