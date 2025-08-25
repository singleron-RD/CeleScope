from celescope.tools.starsolo import (
    Starsolo as tools_Starsolo,
    get_opts_starsolo as tools_opts,
    Mapping,
    Cells,
    Demultiplexing,
)


class Starsolo(tools_Starsolo):
    def __init__(self, args):
        super().__init__(args)
        self.extra_starsolo_args += " --soloStrand Reverse "


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
