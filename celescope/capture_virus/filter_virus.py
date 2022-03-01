
from celescope.tools.capture.filter import Filter, get_opts_filter


def get_opts_filter_virus(parser, sub_program):

    get_opts_filter(parser, sub_program)

def filter_virus(args):

    with Filter_virus(args, display_title='Filtering') as runner:
        runner.run()


class Filter_virus(Filter):
    """

    ## Features
    - Filter background virus UMIs base on a UMI threshold.
    There are three methods to determine the UMI threshold:
        - 'auto' : Using a method similar to cell calling method.
        - 'otsu' : UMI counts are first log 2 transformed and then the threshold is determined by [Otsu's method](https://en.wikipedia.org/wiki/Otsu%27s_method)
        - 'hard' : Using User provided UMI threshold.

    ## Output
    - {sample}_filtered_UMI_tsne.csv : Filtered virus UMI counts of each cell barcode.
    """
