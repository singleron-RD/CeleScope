import pandas as pd

from celescope.flv_trust4.summarize import Summarize
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.tools.plotly_plot import Bar_plot
from celescope.flv_trust4.__init__ import CHAIN, PAIRED_CHAIN
from collections import OrderedDict


def gen_vj_annotation_metrics(df, seqtype):
    """
    Generate vdj Annotation Metrics from contig annotations file.
    """

    def get_vj_spanning_pair():
        """
        Get Productive V-J Spanning_Pair metric from annotation file
        Return productive chain pair number. eg: TRA/TRB or IGH/IGL, IGH/IGK.
        """
        df_productive = df[df["productive"]]

        if seqtype == "BCR":
            df_chain_heavy = df_productive[(df_productive["chain"] == "IGH")]
            df_chain_light = df_productive[
                (df_productive["chain"] == "IGL") | (df_productive["chain"] == "IGK")
            ]
        else:
            df_chain_heavy = df_productive[df_productive["chain"] == "TRA"]
            df_chain_light = df_productive[df_productive["chain"] == "TRB"]

        for _df in [df_chain_heavy, df_chain_light]:
            _df.drop_duplicates(["barcode"], inplace=True)

        VJ_Spanning_Pair_Cells = pd.merge(
            df_chain_heavy, df_chain_light, on="barcode", how="inner"
        )

        return VJ_Spanning_Pair_Cells.shape[0]

    metric_dict = OrderedDict()
    chains, chain_pairs = CHAIN[seqtype], PAIRED_CHAIN[seqtype]

    metric_dict["Cells match with scRNA-seq analysis"] = len(set(df.barcode))
    metric_dict["Cells With Productive V-J Spanning Pair"] = get_vj_spanning_pair()

    for pair in chain_pairs:
        chain1, chain2 = pair.split("_")[0], pair.split("_")[1]
        cbs1 = set(df[(df["productive"]) & (df["chain"] == chain1)].barcode)
        cbs2 = set(df[(df["productive"]) & (df["chain"] == chain2)].barcode)
        metric_dict[f"Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair"] = (
            len(cbs1.intersection(cbs2))
        )

    for chain in chains:
        metric_dict[f"Cells With {chain} Contig"] = len(
            set(df[df["chain"] == chain].barcode)
        )
        metric_dict[f"Cells With CDR3-annotated {chain} Contig"] = len(
            set(df[(df["chain"] == chain) & (df["cdr3"] != "None")].barcode)
        )
        metric_dict[f"Cells With V-J Spanning {chain} Contig"] = len(
            set(df[(df["full_length"]) & (df["chain"] == chain)].barcode)
        )
        metric_dict[f"Cells With Productive {chain} Contig"] = len(
            set(df[(df["productive"]) & (df["chain"] == chain)].barcode)
        )

    return metric_dict


def gen_clonotypes_table(df, out_clonotypes, seqtype):
    """
    Generate clonotypes.csv file
    """
    df = df[df["productive"]]
    df["chain_cdr3aa"] = df[["chain", "cdr3"]].apply(":".join, axis=1)
    df = df.rename(
        columns={"chain_cdr3aa": "cdr3s_aa", "raw_clonotype_id": "clonotype_id"}
    )
    df = df.dropna(subset=["clonotype_id"])
    df = df.sort_values(
        "clonotype_id", key=lambda x: x.str.lstrip("clonotype").astype(int)
    )

    sort_method = {"TCR": True, "BCR": False}
    cdr3_aa_dict = df.groupby("clonotype_id")["cdr3s_aa"].apply(set).to_dict()
    cdr3_aa_dict = {
        key: ";".join(sorted(list(value), reverse=sort_method[seqtype]))
        for key, value in cdr3_aa_dict.items()
    }
    count_dict = df.groupby("clonotype_id")["barcode"].nunique().to_dict()

    df["frequency"] = df["clonotype_id"].apply(lambda x: count_dict[x])
    df["cdr3s_aa"] = df["clonotype_id"].apply(lambda x: cdr3_aa_dict[x])
    df = df.drop_duplicates("clonotype_id")
    df = df[["clonotype_id", "cdr3s_aa", "frequency"]]

    sum_frequency = df["frequency"].sum()
    df["proportion"] = df["frequency"].apply(lambda x: x / sum_frequency)

    df.to_csv(out_clonotypes, sep=",", index=False)


class Annotation(Step):
    """
    ## Features

    - Output TSNE-plot of Assembled T/B Cells.

    ## Output
    - `05.annotation/{sample}_mapping.pdf` TSNE-plot of Assembled Cells.

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.match_dir = args.match_dir
        self.chains, self.paired_groups = Summarize._parse_seqtype(self.seqtype)
        self.contig_file = f"{args.summarize_out}/{self.sample}_filtered_contig.csv"
        self.clonotype_file = f"{args.summarize_out}/clonotypes.csv"

    def parse_clonotype(self):
        """Generate clonotypes table in html."""
        df_clonotypes = pd.read_csv(self.clonotype_file, sep=",")
        df_clonotypes["ClonotypeID"] = df_clonotypes["clonotype_id"].apply(
            lambda x: x.strip("clonetype")
        )
        df_clonotypes["Frequency"] = df_clonotypes["frequency"]
        df_clonotypes["Proportion"] = df_clonotypes["proportion"].apply(
            lambda x: f"{round(x*100, 2)}%"
        )
        df_clonotypes["CDR3_aa"] = df_clonotypes["cdr3s_aa"].apply(
            lambda x: x.replace(";", "<br>")
        )

        return df_clonotypes

    @staticmethod
    @utils.add_log
    def get_vdj_metric(df, chains, pairs):
        """
        Add vdj metrics in html.
        """
        metric_result = []
        fl_pro_pair_df = pd.DataFrame(df[df["productive"]].barcode.value_counts())
        fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df["barcode"] >= 2]
        cell_nums = len(set(df["barcode"]))

        metric_result.append(
            {
                "name": "Cells With Productive V-J Spanning Pair",
                "value": fl_pro_pair_df.shape[0],
                "total": cell_nums,
            }
        )

        for pair in pairs:
            chain1, chain2 = pair.split("_")[0], pair.split("_")[1]
            cbs1 = set(
                df[
                    (df["full_length"]) & (df["productive"]) & (df["chain"] == chain1)
                ].barcode
            )
            cbs2 = set(
                df[
                    (df["full_length"]) & (df["productive"]) & (df["chain"] == chain2)
                ].barcode
            )
            paired_cbs = len(cbs1.intersection(cbs2))

            metric_result.append(
                {
                    "name": f"Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair",
                    "value": paired_cbs,
                    "total": cell_nums,
                    "help_info": "Fraction of cell-associated barcodes with one productive contig for each chain of the receptor pair.A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region",
                }
            )

        for chain in chains:
            metric_result.append(
                {
                    "name": f"Cells With {chain} Contig",
                    "value": len(set(df[df["chain"] == chain].barcode)),
                    "total": cell_nums,
                    "help_info": f"Fraction of cell-associated barcodes with at least one {chain} contig annotated as a full or partial V(D)J gene",
                }
            )
            metric_result.append(
                {
                    "name": f"Cells With CDR3-annotated {chain} Contig",
                    "value": len(
                        set(
                            df[
                                (df["chain"] == chain) & (df["cdr3"] is not None)
                            ].barcode
                        )
                    ),
                    "total": cell_nums,
                }
            )
            metric_result.append(
                {
                    "name": f"Cells With V-J Spanning {chain} Contig",
                    "value": len(
                        set(df[(df["full_length"]) & (df["chain"] == chain)].barcode)
                    ),
                    "total": cell_nums,
                    "help_info": f"Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for {chain}",
                }
            )
            metric_result.append(
                {
                    "name": f"Cells With Productive {chain} Contig",
                    "value": len(
                        set(
                            df[
                                (df["full_length"])
                                & (df["productive"])
                                & (df["chain"] == chain)
                            ].barcode
                        )
                    ),
                    "total": cell_nums,
                    "help_info": "Fraction of cell-associated barcodes with productive IGL chain. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region",
                }
            )

        return metric_result

    @utils.add_log
    def annotation_process(self):
        """Add metrics, clonotypes table, and bar-plot of clonotypes distribution in html."""
        df = pd.read_csv(self.contig_file)

        annotation_summary = Annotation.get_vdj_metric(
            df, self.chains, self.paired_groups
        )
        for anno in annotation_summary:
            self.add_metric(
                anno["name"], anno["value"], anno.get("total"), anno.get("help_info")
            )

        df_clonotypes = self.parse_clonotype()
        title = "Clonetypes"
        table_dict = self.get_table_dict(
            title=title,
            table_id="clonetypes",
            df_table=df_clonotypes[
                ["ClonotypeID", "CDR3_aa", "Frequency", "Proportion"]
            ],
        )
        self.add_data(table_dict=table_dict)

        df_clonotypes["ClonotypeID"] = df_clonotypes["ClonotypeID"].astype("int")
        df_clonotypes.sort_values(by=["ClonotypeID"], inplace=True)
        Barplot = Bar_plot(df_bar=df_clonotypes).get_plotly_div()
        self.add_data(Barplot=Barplot)

    def run(self):
        self.annotation_process()


@utils.add_log
def annotation(args):
    with Annotation(args, display_title="V(D)J Annotation") as runner:
        runner.run()


def get_opts_annotation(parser, sub_program):
    parser.add_argument(
        "--seqtype", help="TCR or BCR", choices=["TCR", "BCR"], required=True
    )
    if sub_program:
        s_common(parser)
        parser.add_argument(
            "--match_dir", help="scRNA-seq match directory", required=True
        )
        parser.add_argument(
            "--summarize_out", help="summarize output directory", required=True
        )
    return parser
