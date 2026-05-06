from celescope.tools import analysis_wrapper
from celescope.tools import utils


@utils.add_log
def analysis(args):
    with analysis_wrapper.Scanpy_wrapper(
        args, display_title="Analysis"
    ) as scanpy_wrapper:
        scanpy_wrapper.run()

    if args.celltypist_model:
        with analysis_wrapper.Celltypist_wrapper(
            args, display_title="Celltypist"
        ) as celltypist_wrapper:
            celltypist_wrapper.run()


def get_opts_analysis(parser, sub_program):
    analysis_wrapper.get_opts_analysis(parser, sub_program)
