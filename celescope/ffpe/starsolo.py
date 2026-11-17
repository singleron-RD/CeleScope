from celescope.tools.starsolo import (
    Starsolo as tools_Starsolo,
    get_opts_starsolo as tools_opts,
    Mapping as tools_Mapping,
    Cells,
    Demultiplexing,
)
from celescope.tools import utils


class Starsolo(tools_Starsolo):
    def __init__(self, args):
        super().__init__(args)
        self.extra_starsolo_args += " --soloStrand Reverse "


class Mapping(tools_Mapping):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    @utils.add_log
    def add_star_metrics(self):
        log_dict = utils.parse_star_log_file(
            f"{self.outdir}/{self.sample}_Log.final.out"
        )
        total_reads = log_dict["Number of input reads"]

        self.add_metric(
            name="Average input read length",
            value=log_dict["Average input read length"],
            help_info="Average length of the input reads after trimming.",
        )
        self.add_metric(
            name="Reads unmapped: too many mismatches",
            value=log_dict["Number of reads unmapped: too many mismatches"],
            total=total_reads,
            help_info="The number of mismatches or the mismatch rate in the best alignment exceeds the allowed threshold. This may be related to sequencing errors, poor sequencing quality, or differences between the sample and the reference genome.",
        )
        self.add_metric(
            name="Reads unmapped: too short",
            value=log_dict["Number of reads unmapped: too short"],
            total=total_reads,
            help_info="The read does not reach the minimum aligned length or the minimum alignment score required by STAR. This may be because the remaining effective sequence after poly-A trimming is too short, or because only a small fraction of the read can match the reference genome.",
        )
        self.add_metric(
            name="Reads unmapped: other",
            value=log_dict["Number of reads unmapped: other"],
            total=total_reads,
            help_info="Reads unmapped for reasons other than the two categories above, usually because STAR cannot find an acceptable alignment. Possible reasons include poor sequence quality, contamination, or sequences that are absent from the reference genome.",
        )

    @utils.add_log
    def run(self):
        metrics = self.parse_cellReadsStats()
        result = self.add_metrics_to_report(metrics)
        self.add_star_metrics()
        return result


def get_opts_starsolo(parser, sub_program=True):
    tools_opts(parser, sub_program)
    parser.set_defaults(
        outFilterMatchNmin=30,
        soloFeatures="GeneFull_Ex50pAS Gene SJ",
    )
    return parser


def starsolo(args):
    with Starsolo(args) as runner:
        q30_cb, q30_umi, chemistry = runner.run()

    with Mapping(args) as runner:
        valid_reads, corrected = runner.run()

    with Cells(args) as runner:
        n_reads, q30_RNA = runner.run(chemistry, valid_reads)

    with Demultiplexing(args) as runner:
        runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)
