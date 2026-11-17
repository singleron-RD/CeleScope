from celescope.tools.tag.mapping_tag import (
    Mapping_tag as Mt,
    get_opts_mapping_tag as opts,
    mapping_tag as mt,
)


# L23C15, not L25C15 like Clindex
def get_opts_mapping_tag(parser, sub_program):
    opts(parser, sub_program)
    parser.set_defaults(fq_pattern="L23C15")


class Mapping_tag(Mt):
    pass


def mapping_tag(args):
    mt(args)
