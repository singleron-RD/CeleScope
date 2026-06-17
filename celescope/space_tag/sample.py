from celescope.tools.sample import (
    Sample as Sa,
    get_opts_sample as opts,
    sample as sa,
)


class Sample(Sa):
    pass


def sample(args):
    sa(args)


def get_opts_sample(parser, sub_program):
    opts(parser, sub_program)
