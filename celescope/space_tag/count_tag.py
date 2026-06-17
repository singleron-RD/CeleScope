import numpy as np

from celescope.tools.tag.count_tag import (
    Count_tag as Ct,
    get_opts_count_tag as opts,
)


class Count_tag(Ct):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)

    def add_seq_metrics(self):
        self.add_metric(
            name="Number of Positive Tissue Spots",
            value=self.n_match_barcode,
            help_info="Number of in-tissue Spots with as least one tag UMI",
        )

        self.add_metric(
            name="Mapped Reads in Spots",
            value=self.mapped_read_in_cell,
            total=self.mapped_read,
            help_info="Mapped reads with in-tissue Spots",
        )

        UMIs = self.df_UMI_cell.apply(sum, axis=1)
        umi_median = round(np.median(UMIs), 2)
        umi_mean = round(np.mean(UMIs), 2)
        self.add_metric(
            name="Median UMI per Spots",
            value=umi_median,
            help_info="Median UMI per in-tissue Spot",
        )
        self.add_metric(
            name="Mean UMI per Spots",
            value=umi_mean,
            help_info="Mean UMI per in-tissue Spot",
        )


def count_tag(args):
    with Count_tag(args) as runner:
        runner.run()


def get_opts_count_tag(parser, sub_program):
    opts(parser, sub_program)
