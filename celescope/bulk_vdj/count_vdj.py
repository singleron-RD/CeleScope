import pandas as pd

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj.__init__ import CHAINS


class Count_vdj(Step):
    """
    ## Features  
    - Summarize clonetypes infomation.

    ## Output
    - `{sample}_filtered_annotations.csv` Annotations for each CDR3.

    - `{sample}_clonetypes.csv` The count and percentage of each CDR3 of reads.

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.chains = CHAINS[args.type]

        # IN
        self.productive_file = args.productive_file

        # OUT
        self.annotation = f"{self.out_prefix}_filtered_annotations.csv"
        self.clonetypes = f"{self.out_prefix}_clonetypes.csv"

    
    @utils.add_log
    def gen_clonotypes(self):
        """
        Output clonotypes.csv file.

        CDR3_ID	aaSeqCDR3	Frequency	Proportion
        1	TRB:CAILDGRYEQYF	1363439	12.93%
        2	TRB:CASRTTDWRDYGYTF	355253	3.37%
        3	TRA:CAVAHHTGTASKLTF	336672	3.19%
        """
        productive_file = pd.read_csv(self.productive_file, sep='\t')
        productive_file = productive_file[productive_file["chain"].isin(self.chains)]

        groupby_elements = [
            'chain',
            'aaSeqCDR3',
        ]
        clonetypes = productive_file.groupby(groupby_elements, as_index=False).agg({"readID": "count"})
        clonetypes["aaSeqCDR3"] = clonetypes.loc[:, ["chain", "aaSeqCDR3"]].apply(':'.join, axis=1)
        clonetypes = clonetypes.sort_values("readID", ascending=False)
        clonetypes = clonetypes.reset_index()
        clonetypes["CDR3_ID"] = pd.Series(clonetypes.index) + 1
        clonetypes = clonetypes.rename(columns={"readID": "Frequency"})

        sum_frequency = sum(clonetypes["Frequency"])
        clonetypes["Proportion"] = clonetypes["Frequency"].apply(lambda x : x / sum_frequency * 100)
        clonetypes["Proportion"] = clonetypes["Proportion"].apply(lambda x: str(round(x, 2)) + '%' )
        
        clonetypes = clonetypes[["CDR3_ID", "aaSeqCDR3", "Frequency", "Proportion"]]
        clonetypes.to_csv(self.clonetypes, sep=',', index=False)

        top_500_clonetypes = clonetypes.head(500)
        table_dict = self.get_table_dict(
            title="Top500_CDR3",
            table_id='clonetypes',
            df_table=top_500_clonetypes
        )
        self.add_data(table_dict=table_dict)

        CDR3_dict = dict(zip(clonetypes["aaSeqCDR3"], clonetypes["CDR3_ID"]))
        return CDR3_dict

    @utils.add_log
    def gen_annotation(self, CDR3_dict):
        """
        Output annotation file

        readID,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,raw_cdr3_id
        readId_1,TRB,TRBV11-2,None,TRBJ1-1,None,True,True,CASSSDRGWSEAFF,TGTGCCAGCAGCTCCGACAGGGGTTGGTCTGAAGCTTTCTTT,10
        readId_5,TRB,TRBV18,None,TRBJ1-4,None,True,True,CASSPGENNEKLFF,TGTGCCAGCTCCCCTGGGGAGAATAATGAAAAACTGTTTTTT,11
        readId_7,TRB,TRBV7-9,None,TRBJ2-7,None,True,True,CAILDGRYEQYF,TGTGCCATCCTAGACGGGCGCTACGAGCAGTACTTC,1
        """
        in_file = pd.read_csv(self.productive_file, sep='\t')
        in_file = in_file[in_file["chain"].isin(self.chains)]
        in_file["chain_cdr3"] = in_file.loc[:, ["chain", "aaSeqCDR3"]].apply(':'.join, axis=1)
        out_file = pd.DataFrame({
                        'readID':in_file.readID,
                        'chain':in_file.chain,
                        'v_gene':in_file.bestVGene,
                        'd_gene':'None',
                        'j_gene':in_file.bestJGene,
                        'c_gene':'None',
                        'full_length':True,
                        'productive':True,
                        'cdr3':in_file.aaSeqCDR3,
                        'cdr3_nt':in_file.nSeqCDR3,
                        'raw_cdr3_id':in_file.chain_cdr3.apply(lambda x: CDR3_dict[x])
                        })
        out_file.to_csv(self.annotation, sep=',', index=False)

    def run(self):
        CDR3_dict = self.gen_clonotypes()
        self.gen_annotation(CDR3_dict)


def count_vdj(args):
    with Count_vdj(args, display_title="Count") as runner:
        runner.run()


def get_opts_count_vdj(parser, sub_program):
    parser.add_argument(
        "--type", help="Required. `TCR` or `BCR`. ", required=True)
    if sub_program:
        parser.add_argument("--productive_file", help="Required. File from step mapping_vdj.", required=True)
        parser = s_common(parser)
