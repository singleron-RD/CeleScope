from celescope.tools.starsolo import (
    Starsolo,
    Mapping,
    Cells,
    Demultiplexing,
    get_opts_starsolo as opts_super,
)


def starsolo(args):
    args.STAR_param += " --outSAMunmapped Within "
    with Starsolo(args) as runner:
        q30_cb, q30_umi, chemistry = runner.run()
    with Mapping(args, display_title="Mapping to Host Genome") as runner:
        valid_reads, corrected = runner.run()
    with Cells(args) as runner:
        n_reads, q30_RNA = runner.run(chemistry, valid_reads)
    with Demultiplexing(args) as runner:
        runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)


def get_opts_starsolo(parser, sub_program):
    opts_super(parser, sub_program)
