import itertools
from celescope.tools.utils import add_log
from celescope.tools.Step import Step, s_common


@add_log
def target_metrics(args):
    step_name = "target_metrics"
    step = Step(args, step_name)

    

    step.clean_up()


def get_opts_target_metrics(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument("--bam", help='featureCounts bam', required=True)
    parser.add_argument("--gene_list", help='gene_list', required=True)

