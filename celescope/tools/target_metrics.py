import itertools
import pysam

from celescope.tools.utils import add_log, read_one_col
from celescope.tools.Step import Step, s_common


@add_log
def target_metrics(args):
    step_name = "target_metrics"
    step = Step(args, step_name)

    gene_list, n_gene = read_one_col(args.gene_list)
    step.add_metric(
        name="Number of Target Genes",
        value=n_gene,
    )

    with pysam.AlignmentFile(args.bam, "rb") as reader:
        for record in reader:
            gene_name = record.get_tag('GN')
            barcode = record.get_tag('CB')
            UMI = record.get_tag('UB')
            

    step._clean_up()


def get_opts_target_metrics(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument("--bam", help='featureCounts bam', required=True)
    parser.add_argument("--gene_list", help='gene_list', required=True)

