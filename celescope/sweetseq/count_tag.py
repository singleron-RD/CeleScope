from celescope.tools.tag.count_tag import (
    Count_tag as Ct,
    count_tag as ct,
    get_opts_count_tag as opts,
)


class Count_tag(Ct):
    pass


def count_tag(args):
    ct(args)


def get_opts_count_tag(parser, sub_program):
    opts(parser, sub_program)
