from celescope.tools.starsolo import (
    Starsolo as tools_starsolo,
    Mapping as tools_mapping,
    Cells as tools_cells,
    Demultiplexing,
    get_opts_starsolo as opts_super,
    COUNTS_FILE_NAME,
)


class Starsolo(tools_starsolo):
    def __init__(self, args):
        super().__init__(args)
        self.outs = []


class Mapping(tools_mapping):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        self.outs = []
        solo_dir = f"{self.outdir}/{self.sample}_Solo.out/{args.report_soloFeature}"
        self.filtered_matrix = f"{solo_dir}/filtered"


class Cells(tools_cells):
    def __init__(self, args):
        super().__init__(args)
        solo_dir = solo_dir = (
            f"{self.outdir}/{self.sample}_Solo.out/{args.report_soloFeature}"
        )
        self.counts_file = f"{solo_dir}/{COUNTS_FILE_NAME }"


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
