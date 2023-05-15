import pandas as pd

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj.__init__ import CHAINS
from celescope.tools.plotly_plot import Bar_plot

def correct_cdr3_nt(umi_dict, percent=0.1):
    """
    Correct umi_dict in place.
    Args:
        umi_dict: {cdr3_nt: umi_count}
        percent: if hamming_distance(low_umi_cdr3, high_umi_cdr3) == 1 and
            low_count / high_count < percent, merge low to high.
    Returns:
        correct_dict: dict {low_umi_cdr3: high_umi_cdr3}
    """
    correct_dict = dict()
    
    umi_arr = sorted(
        umi_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    
    while True:
        # break when only highest in umi_arr
        if len(umi_arr) == 1:
            break
        umi_low = umi_arr.pop()
        low_seq = umi_low[0]
        low_count = umi_low[1]

        for umi_kv in umi_arr:
            high_seq = umi_kv[0]
            high_count = umi_kv[1]
            if len(low_seq) != len(high_seq):
                break
            if float(low_count / high_count) > percent:
                break
            if utils.hamming_distance(low_seq, high_seq) == 1:
                correct_dict[low_seq] = high_seq
                n_low = umi_dict[low_seq]
                # merge
                umi_dict[high_seq] += n_low
                del (umi_dict[low_seq])
                break
            
    return correct_dict


class Count_vdj(Step):
    """
    ## Features  
    - Summarize clonetypes infomation.
    ## Output
    - `{sample}_corrected_productive.tsv` All productive chains info after cdr3_nt correcting.
    - `{sample}_filtered_annotations.csv` Annotations for each CDR3.
    - `{sample}_clonetypes.csv` The count and percentage of each CDR3_aa of umis.
    - `{sample}_clonetypes_by_nt.csv` The count and percentage of each CDR3_nt of umis.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.chains = CHAINS[args.type]

        # IN
        self.productive_file = args.productive_file

        # OUT
        self.corrected_productive_file = f"{self.out_prefix}_corrected_productive.tsv"
        self.annotation = f"{self.out_prefix}_filtered_annotations.csv"
        self.clonetypes_aa = f"{self.out_prefix}_clonetypes.csv"
        self.clonetypes_nt = f"{self.out_prefix}_clonetypes_by_nt.csv"
        
    @utils.add_log
    def correct_cdr3_nt(self):
        """
        Correct cdr3_nt sequence.
        """
        self.productive_file = pd.read_csv(self.productive_file, sep='\t')
        self.productive_file = self.productive_file[self.productive_file["chain"].isin(self.chains)]
        self.productive_file["umi"] = self.productive_file["sequence_id"].apply(lambda x: x.split('_')[1])

        for chain in self.chains:
            df_tmp = self.productive_file[self.productive_file["chain"] == chain]
            groupby_elements = ["chain","nSeqCDR3"]
            clonetypes = df_tmp.groupby(groupby_elements, as_index=False).agg({"umi": "count"})
            clonetypes = clonetypes.sort_values("umi", ascending=False)
            umi_dict = dict(zip(list(clonetypes.nSeqCDR3),list(clonetypes.umi)))
            correct_dict = correct_cdr3_nt(umi_dict)
    
            for low_nt, high_nt in correct_dict.items():
                high_aa = df_tmp.loc[df_tmp.nSeqCDR3 == high_nt].aaSeqCDR3.iloc[0]
                self.productive_file.loc[self.productive_file.nSeqCDR3 == low_nt, 'aaSeqCDR3'] = high_aa
                self.productive_file.loc[self.productive_file.nSeqCDR3 == low_nt, 'nSeqCDR3'] = high_nt
                
        self.productive_file.to_csv(self.corrected_productive_file, sep="\t", index=False)
    
    @utils.add_log
    def gen_clonotypes(self):
        """
        Output clonotypes.csv file.
        
        ClonotypeID,aaSeqCDR3,Frequency,Proportion
        1,TRB:CASSPIRDGVYGYTF,135,0.58%
        2,TRB:CASSEATTQYF,105,0.45%
        3,TRA:CALSESKTSYDKVIF,103,0.44%
        4,TRA:CAVSSSDRGSTLGRLYF,98,0.42%
        5,TRB:CSVAQGYGDTEAFF,95,0.41%
        
        """
        CDR3_seq = ["nSeqCDR3", "aaSeqCDR3"]
        for seq in CDR3_seq:
            groupby_elements = ["chain", seq]
            #clonetypes = productive_file.groupby(groupby_elements, as_index=False).agg({"sequence_id": "count"})
            clonetypes = self.productive_file.groupby(groupby_elements, as_index=False).agg({"umi": "count"})
            clonetypes[seq] = clonetypes.loc[:, ["chain", seq]].apply(':'.join, axis=1)
            clonetypes = clonetypes.sort_values("umi", ascending=False)
            clonetypes = clonetypes.reset_index()
            clonetypes["ClonotypeID"] = pd.Series(clonetypes.index) + 1
            clonetypes = clonetypes.rename(columns={"umi": "Frequency"})

            sum_frequency = sum(clonetypes["Frequency"])
            clonetypes["Proportion"] = clonetypes["Frequency"].apply(lambda x : x / sum_frequency)
            proportion_list = clonetypes["Proportion"].tolist()
            clonetypes["Proportion"] = clonetypes["Proportion"].apply(lambda x: str(round(x*100, 2)) + '%' )
            clonetypes = clonetypes[["ClonotypeID", seq, "Frequency", "Proportion"]]
            
            if seq == "aaSeqCDR3":
                clonetypes.to_csv(self.clonetypes_aa, sep=',', index=False)
                CDR3_dict = dict(zip(clonetypes[seq], clonetypes["ClonotypeID"]))
                top_500_clonetypes = clonetypes.head(500)
                table_dict = self.get_table_dict(
                    title="Top500",
                    table_id='clonetypes',
                    df_table=top_500_clonetypes
                )
                self.add_data(table_dict=table_dict)
                
                clonetypes["CDR3_aa"] = clonetypes[seq].apply(lambda x: x.replace(';', '<br>'))
                clonetypes["proportion"] = proportion_list
                Barplot = Bar_plot(df_bar=clonetypes).get_plotly_div()
                self.add_data(Barplot=Barplot)
            else:
                clonetypes.to_csv(self.clonetypes_nt, sep=',', index=False)
                
        return CDR3_dict

    @utils.add_log
    def gen_annotation(self, CDR3_dict):
        """
        Output annotation file
        """
        in_file = self.productive_file
        in_file["chain_cdr3"] = in_file.loc[:, ["chain", "aaSeqCDR3"]].apply(':'.join, axis=1)
        umi_dict = in_file.groupby(["chain_cdr3"], as_index=False).agg({"umi":"count"}).set_index("chain_cdr3").to_dict()["umi"]
        in_file["reads"] = in_file["chain_cdr3"].apply(lambda x: umi_dict[x])
        in_file["umi"] = in_file["chain_cdr3"].apply(lambda x: umi_dict[x])

        out_file = pd.DataFrame({
                        "barcode":in_file.barcode,
                        "is_cell":True,
                        "sequence_id":in_file.sequence_id,
                        "high_confidence":True,
                        "chain":in_file.chain,
                        "v_gene":in_file.bestVGene,
                        "d_gene":in_file.bestDGene,
                        "j_gene":in_file.bestJGene,
                        "c_gene":'None',
                        "full_length":True,
                        "productive":True,
                        "cdr3":in_file.aaSeqCDR3,
                        "cdr3_nt":in_file.nSeqCDR3,
                        "reads":in_file.reads,
                        "umis":in_file.umi,
                        "raw_cdr3_id":in_file.chain_cdr3.apply(lambda x: CDR3_dict[x])
                        })
        out_file.fillna("None", inplace=True)
        out_file.to_csv(self.annotation, sep=',', index=False)

    def run(self):
        self.correct_cdr3_nt()
        CDR3_dict = self.gen_clonotypes()
        self.gen_annotation(CDR3_dict)


def count_vdj(args):
    with Count_vdj(args, display_title="Clonotypes") as runner:
        runner.run()


def get_opts_count_vdj(parser, sub_program):
    parser.add_argument(
        "--type", help="Required. `TCR` or `BCR`. ", required=True)
    if sub_program:
        parser.add_argument("--productive_file", help="Required. File from step mapping_vdj.", required=True)
        parser = s_common(parser)