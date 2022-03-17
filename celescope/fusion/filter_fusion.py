
from celescope.tools.capture.filter import Filter, get_opts_filter


def get_opts_filter_fusion(parser, sub_program):

    get_opts_filter(parser, sub_program)

def filter_fusion(args):

    with Filter_fusion(args) as runner:
        runner.run()


class Filter_fusion(Filter):
    """    
    ## Features
    - Correct single-base errors in UMIs due to sequencing, amplification, etc.
    - Filter background UMIs base on a UMI threshold.
    There are three methods to determine the UMI threshold:
        - 'auto' : Using a method similar to cell calling method.
        - 'otsu' : UMI counts are first log 2 transformed and then the threshold is determined by [Otsu's method](https://en.wikipedia.org/wiki/Otsu%27s_method)
        - 'hard' : Using User provided UMI threshold.

    ## Output
    - `{sample}_corrected_read_count.json` Read counts after UMI correction.
    - `{sample}_filtered_read_count.json` Filtered read counts.
    - `{sample}_filtered_UMI.csv` Filtered UMI counts.

    """