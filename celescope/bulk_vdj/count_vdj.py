import pandas as pd
import numpy as np
import numbers
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj.__init__ import CHAINS


def format_value(value, total):
    if not isinstance(value, numbers.Number):
        return value
    display = str(format(value, ','))
    if total:
        fraction = round(value / total * 100, 2)
        display += f'({fraction}%)'
    return display


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


def simpson_di(data):

    """ Given a hash { 'species': count } , returns the Simpson Diversity Index
    
    >>> simpson_di({'a': 10, 'b': 20, 'c': 30,})
    0.3888888888888889
    """

    def p(n, N):
        """ Relative abundance """
        if n == 0:
            return 0
        else:
            return float(n)/N

    N = sum(data.values())
    
    return sum(p(n, N)**2 for n in data.values() if n != 0)


def inverse_simpson_di(data):
    """ Given a hash { 'species': count } , returns the inverse Simpson Diversity Index
    
    >>> inverse_simpson_di({'a': 10, 'b': 20, 'c': 30,})
    2.571428571428571
    """
    return float(1)/simpson_di(data)


class Count_vdj(Step):
    """
    ## Features  
    - Summarize clonetypes infomation.
    ## Output
    - `{sample}_mapping_metrics.tsv` Mapping metrics of each index.
    - `{sample}_corrected_productive.tsv` Productive chains after cdr3_nt correcting.
    - `{sample}_filtered_annotations.csv` Annotation for each CDR3.
    - `{sample}_clonetypes.csv` The count and percentage of each CDR3_aa of umis.
    - `{sample}_clonetypes_by_nt.csv` The count and percentage of each CDR3_nt of umis.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.chains = CHAINS[args.type]

        # IN
        self.airr_file = args.airr_file
        self.productive_file = args.productive_file

        # OUT
        self.mapping_result_file = f"{self.out_prefix}_mapping_metrics.tsv"
        self.corrected_productive_file = f"{self.out_prefix}_corrected_productive.tsv"
        self.annotation = f"{self.out_prefix}_filtered_annotations.csv"
        self.clonetypes_aa = f"{self.out_prefix}_clonetypes.csv"
        self.clonetypes_nt = f"{self.out_prefix}_clonetypes_nt.csv"
    
    @utils.add_log
    def gen_mapping_metrics(self):
        """mapping metrics for each index.
        """
        df = pd.read_csv(self.airr_file, sep='\t')
        df.fillna("", inplace=True)
        df["index"] = df["sequence_id"].apply(lambda x: x.split('_')[0])
        
        df_metrics = pd.DataFrame(columns=[
            "Index", "UMI_Counts", "UMIs_Mapped_To_Any_VDJ_Gene", "UMIs_with_CDR3",
            "UMIs_with_Correct_CDR3", "UMIs_Mapped_Confidently_To_VJ_Gene",
        ])
        for chain in self.chains:
            df_metrics[f"UMIs_Mapped_To_{chain}"] = ''
            
        for index, data in df.groupby("index"):
            umi_count = data.shape[0]
            
            data = data[(data["v_call"]!="") | ((data["d_call"]!="")) | ((data["j_call"]!=""))]
            UMIs_Mapped_To_Any_VDJ_Gene = data.shape[0]

            df_cdr3 = data[(data["cdr3_aa"]!="") & (data["junction_aa"]!="")]
            UMIs_with_CDR3 = df_cdr3.shape[0]
            
            df_correct_cdr3 = df_cdr3[~(df_cdr3["cdr3_aa"].str.contains(r"\*")) & ~(df_cdr3["cdr3_aa"].str.contains("X"))]
            UMIs_with_Correct_CDR3 = df_correct_cdr3.shape[0]
            
            df_confident = df_correct_cdr3[df_correct_cdr3["productive"]=="T"]
            UMIs_Mapped_Confidently_To_VJ_Gene = df_confident.shape[0]
            
            UMIs_Mapped_To_chains = []
            for chain in self.chains:
                df_chain = df_confident[df_confident.locus == chain]
                UMIs_Mapped_To_chains.append(df_chain.shape[0])
            
            metrics_list = [
                index, umi_count, UMIs_Mapped_To_Any_VDJ_Gene, UMIs_with_CDR3,
                UMIs_with_Correct_CDR3, UMIs_Mapped_Confidently_To_VJ_Gene,
            ]
            metrics_list += UMIs_Mapped_To_chains
            metrics_list = [format_value(value, umi_count) for value in metrics_list]
            
            df_metrics.loc[len(df_metrics.index)] = metrics_list
            df_metrics.to_csv(self.mapping_result_file, sep="\t", index=False)
        
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
    def gen_clonotype_file(self):
        """Output clonotypes.csv file.  
        
        ClonotypeID	aaSeqCDR3	Frequency	Proportion	Index	diversity
        1	TRB:CASGDRLGGSQNTLYF	105	2.34%	ACTTGT	463.14
        2	TRA:CAASANSGTYQRF	36	0.80%	ACTTGT	463.14
        3	TRA:CALSRNTGYQNFYF	30	0.67%	ACTTGT	463.14
        4	TRA:CAASADYGNEKITF	30	0.67%	ACTTGT	463.14
        5	TRA:CAASVDNYAQGLTF	27	0.60%	ACTTGT	463.14
        """
        index_set = set(self.productive_file.barcode)
        cdr3_types = ["nSeqCDR3", "aaSeqCDR3"]
        final_df_aa = pd.DataFrame()
        final_df_nt = pd.DataFrame()
        
        def format_files(index, cdr3_type):
            """format clonotypes.csv and clonotypes_nt.csv from productive file."""
            nonlocal final_df_aa, final_df_nt
            
            groupby_elements = ["chain", cdr3_type]
            df_clonetypes = self.productive_file[self.productive_file["barcode"]==index]
            df_clonetypes = df_clonetypes.groupby(groupby_elements, as_index=False).agg({"umi": "count"})
            df_clonetypes[cdr3_type] = df_clonetypes.loc[:, ["chain", cdr3_type]].apply(':'.join, axis=1)
            df_clonetypes = df_clonetypes.sort_values("umi", ascending=False).reset_index()
            df_clonetypes["ClonotypeID"] = pd.Series(df_clonetypes.index) + 1
            df_clonetypes = df_clonetypes.rename(columns={"umi": "Frequency"})
    
            sum_frequency = sum(df_clonetypes["Frequency"])
            df_clonetypes["Proportion"] = df_clonetypes["Frequency"].apply(lambda x : x / sum_frequency)
            df_clonetypes["Proportion"] = df_clonetypes["Proportion"].apply(lambda x: str(round(x*100, 2)) + '%' )
            df_clonetypes = df_clonetypes[["ClonotypeID", cdr3_type, "Frequency", "Proportion"]]
            df_clonetypes["Index"] = index

            # calculate diversity
            data = dict(zip(df_clonetypes[cdr3_type], df_clonetypes["Frequency"]))
            clonotype_diversity = round(inverse_simpson_di(data), 2)
            df_clonetypes["Diversity"] = clonotype_diversity
            # df_clonetypes.loc[df_clonetypes.Index == index, "Diversity"] = clonotype_diversity
            if cdr3_type == "aaSeqCDR3":
                final_df_aa = pd.concat([final_df_aa, df_clonetypes])
            else:
                final_df_nt = pd.concat([final_df_nt, df_clonetypes])

        for index in index_set:
            for cdr3_type in cdr3_types:
                format_files(index, cdr3_type)

        final_df_aa.to_csv(self.clonetypes_aa, sep=',', index=False)
        final_df_nt.to_csv(self.clonetypes_nt, sep=',', index=False)
        
        df_table = final_df_aa.groupby('Index').head(100)
        df_table = df_table[["Index", "ClonotypeID", "aaSeqCDR3", "Frequency", "Proportion", "Diversity"]]
        
        mean_diversity = round(np.mean(df_table.Diversity), 2)
        self.add_metric(
            name="Mean Diversity",
            value=mean_diversity,
            help_info="Mean diversity of all indexes"
        )
        
        table_dict = self.get_table_dict(
            title="Clonotypes by Index",
            table_id='bulk_vdj',
            df_table=df_table
        )
        self.add_data(table_dict=table_dict)


    @utils.add_log
    def gen_annotation_file(self):
        """
        Output annotation file
        """
        in_file = self.productive_file.copy()
        
        groupby_elements = ["barcode", "chain", "aaSeqCDR3"]
        in_file['umis'] = in_file.groupby(groupby_elements)['umi'].transform('count')

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
                        "reads":in_file.umis,
                        "umis":in_file.umis,
                        })
        out_file.fillna("None", inplace=True)
        out_file.to_csv(self.annotation, sep=',', index=False)

    def run(self):
        self.gen_mapping_metrics()
        self.correct_cdr3_nt()
        self.gen_clonotype_file()
        self.gen_annotation_file()


def count_vdj(args):
    with Count_vdj(args, display_title="Clonotypes") as runner:
        runner.run()


def get_opts_count_vdj(parser, sub_program):
    parser.add_argument(
        "--type", help="Required. `TCR` or `BCR`. ", required=True)
    if sub_program:
        parser.add_argument("--productive_file", help="Required. Productive file from mapping_vdj step.", required=True)
        parser.add_argument("--airr_file", help="Required. Airr file from mapping_vdj step.", required=True)
        parser = s_common(parser)