from celescope.tools.barcode import get_opts_barcode
from celescope.tools.cutadapt import get_opts_cutadapt
from celescope.rna.star import get_opts_star

class prep:
    pass

def get_opts_prep(parser, sub_program):
    get_opts_barcode(parser, sub_program=False)
    get_opts_cutadapt(parser, sub_program=False)
    get_opts_star(parser, sub_program=False)