from celescope.tools.capture.analysis import Analysis, get_opts_analysis


def analysis_virus(args):

    with Analysis_virus(args) as runner:
        runner.run()


def get_opts_analysis_virus(parser, sub_program):
    get_opts_analysis(parser, sub_program)


class Analysis_virus(Analysis):
    pass

