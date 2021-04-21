import itertools
import pysam
import collections
import numpy as np

import celescope.tools.utils as utils
from celescope.tools.Step import Step, s_common


@utils.add_log
def target_metrics(args):
    step_name = "target_metrics"
    step = Step(args, step_name)

    gene_list, n_gene = utils.read_one_col(args.gene_list)
    step.add_metric(
        name="Number of Target Genes",
        value=n_gene,
    )

    match_barcode = set(utils.parse_match_dir(args.match_dir)["match_barcode"])

    count_dict = utils.genDict(dim=3, valType=int)
    with pysam.AlignmentFile(args.bam, "rb") as reader:
        for record in reader:
            try:
                gene_name = record.get_tag('GN')
            except KeyError:
                continue
            barcode = record.get_tag('CB')
            UMI = record.get_tag('UB')
            count_dict[barcode][gene_name][UMI] += 1

    UMIs = 0
    enriched_UMIs = 0
    enriched_UMIs_in_cells = 0
    enriched_UMIs_per_cell = []


    for barcode in count_dict:
        barcode_enriched_UMI = 0
        for gene_name in count_dict[barcode]:
            gene_UMI = len(count_dict[barcode][gene_name])
            UMIs += gene_UMI
            if gene_name in gene_list:
                enriched_UMIs += gene_UMI
                if barcode in match_barcode:
                    enriched_UMIs_in_cells += gene_UMI
                    barcode_enriched_UMI += gene_UMI
        if barcode in match_barcode:
            enriched_UMIs_per_cell.append(barcode_enriched_UMI)
    target_metrics.logger.debug(enriched_UMIs_per_cell)

    step.add_metric(
        name="Enriched UMIs",
        value=enriched_UMIs,
        total=UMIs,
    )
    step.add_metric(
        name="Enriched UMIs in Cells",
        value=enriched_UMIs_in_cells,
        total=UMIs,
    )
    step.add_metric(
        name="Median Enriched UMIs per Cell",
        value=np.median(enriched_UMIs_per_cell),
    )

    step.clean_up()


def get_opts_target_metrics(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument("--bam", help='featureCounts bam', required=True)
        parser.add_argument('--match_dir', help='match_dir', required=True)
    parser.add_argument("--gene_list", help='gene_list', required=True)

