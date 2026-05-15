from celescope.tools import analysis_wrapper
from celescope.tools import utils


@utils.add_log
def analysis(args):
    with analysis_wrapper.Analysis(args) as scanpy_wrapper:
        scanpy_wrapper.run()

    if args.celltypist_model:
        with analysis_wrapper.Celltypist(args) as celltypist_wrapper:
            celltypist_wrapper.run()


def get_opts_analysis(parser, sub_program):
    analysis_wrapper.get_opts_analysis(parser, sub_program)
