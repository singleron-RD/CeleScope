from celescope.tools.barcode import (
    Barcode as Ba,
    get_opts_barcode as opts,
    barcode as ba,
)


class Barcode(Ba):
    pass


def barcode(args):
    ba(args)


def get_opts_barcode(parser, sub_program):
    opts(parser, sub_program)
