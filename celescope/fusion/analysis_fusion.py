from celescope.tools.capture.analysis import Analysis, get_opts_analysis


def analysis_fusion(args):

    with Analysis_fusion(args) as runner:
        runner.run()


def get_opts_analysis_fusion(parser, sub_program):
    get_opts_analysis(parser, sub_program)


class Analysis_fusion(Analysis):
    pass