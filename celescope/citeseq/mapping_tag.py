from celescope.tools.tag.mapping_tag import (
    Mapping_tag as Mt,
    get_opts_mapping_tag as opts,
    mapping_tag as mt,
)


class Mapping_tag(Mt):
    pass


def mapping_tag(args):
    mt(args)


def get_opts_mapping_tag(parser, sub_program):
    opts(parser, sub_program)
