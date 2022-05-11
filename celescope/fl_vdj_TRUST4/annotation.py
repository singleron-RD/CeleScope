import pandas as pd
import numpy as np 
import subprocess
import glob
from Bio.Seq import Seq
from celescope.fl_vdj_TRUST4.summarize import Summarize
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.fl_vdj_TRUST4.__init__ import TOOLS_DIR
from celescope.tools.plotly_plot import Bar_plot


@utils.add_log
def run_mapping(rds, contig, sample, outdir, assign):
    """ Generate assembled barcodes results match with scRNA-seq.
    Args:
        rds: rds path of scRNA-seq.
        contig: filtered contig path.
        outdir: output directory.
        assign: auto-assigned info of scRNA-seq.
    """
    cmd = (
        f'Rscript {TOOLS_DIR}/VDJmapping.R '
        f'--rds {rds} '
        f'--VDJ {contig} '
        f'--sample {sample} '
        f'--outdir {outdir} '
        f'--assign_file {assign}'
    )
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def get_vj_annot(df, chains, pairs):
    """ Generate V(D)J annotation for html
    Args:
        chains: TRA, TRB for TCR or IGH, IGK, IGL for BCR.
        pairs: TRA_TRB for TCR or IGH_IGK, IGH_IGL BCR.
    Returns:
        V(D)J annotation for html.
    """
    fl_pro_pair_df = pd.DataFrame(df[df['productive']==True].barcode.value_counts())
    fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']>=2]
    Result = []
    cell_nums = len(set(df['barcode'].tolist()))
    Result.append({
        'name': 'Cells With Productive V-J Spanning Pair',
        'value': fl_pro_pair_df.shape[0],
        'total': cell_nums,
    })
    for p in pairs:
        chain1 = p.split('_')[0]
        chain2 = p.split('_')[1]
        cbs1 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain1)].barcode.tolist())
        cbs2 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain2)].barcode.tolist())
        paired_cbs = len(cbs1.intersection(cbs2))
        Result.append({
            'name': f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
            'value': paired_cbs,
            'total': cell_nums,
            'help_info': "Fraction of cell-associated barcodes with one productive contig for each chain of the receptor pair.A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
        })
    for c in chains:
        Result.append({
            'name': f'Cells With {c} Contig',
            'value': len(set(df[df['chain']==c].barcode.tolist())),
            'total': cell_nums,
            'help_info': f'Fraction of cell-associated barcodes with at least one {c} contig annotated as a full or partial V(D)J gene'
        })
        Result.append({
            'name': f'Cells With CDR3-annotated {c} Contig',
            'value': len(set(df[(df['chain']==c)&(df['productive']==True)].barcode.tolist())),
            'total': cell_nums,
        })
        Result.append({
            'name': f'Cells With V-J Spanning {c} Contig',
            'value': len(set(df[(df['full_length']==True)&(df['chain']==c)].barcode.tolist())),
            'total': cell_nums,
            'help_info': f"Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for {c}"
        })
        Result.append({
            'name': f'Cells With Productive {c} Contig',
            'value': len(set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==c)].barcode.tolist())),
            'total': cell_nums,
            'help_info': "Fraction of cell-associated barcodes with productive IGL chain. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
        })

    return Result


class Annotation(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.match_dir = args.match_dir
        self.chains, self.paired_groups = Summarize._parse_seqtype(self.seqtype)

        try:
            self.rds = glob.glob(f'{self.match_dir}/06.analysis/*.rds')[0]
            self.assign_file = glob.glob(f'{self.match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')[0]
        except IndexError:
            pass
  
        self.contig = glob.glob(f'{self.outdir}/../04.summarize/{self.sample}_filtered_contig.csv')[0]

        if self.seqtype == 'TCR':
            self.Celltype = {'T_cells','NKT_cells','T cells','NK T cells','Tcells'}
            self._name = "Tcells"
        elif self.seqtype == 'BCR':
            self.Celltype = {'Plasma_cells','B_cells','Mature_B_cell','Plasma cells','B cells','Bcells'}
            self._name = "Bcells"
    
    def parse_clonotype(self):
        df_clonotypes=pd.read_csv(f'{self.outdir}/../04.summarize/clonotypes.csv', sep=',')
        df_clonotypes['ClonotypeID'] = df_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        df_clonotypes['Frequency'] = df_clonotypes['frequency']
        df_clonotypes['Proportion'] = df_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        df_clonotypes['CDR3_aa'] = df_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))
        return df_clonotypes

    def parse_contig(self):
        df = pd.read_csv(f'{self.outdir}/../04.summarize/{self.sample}_contig.csv', sep='\t', header=None) 
        df.columns = ['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis']
        df['d_gene'] = df['d_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['c_gene'] = df['c_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['cdr3'] = df['cdr3_nt'].apply(lambda x: 'None' if "*" in str(Seq(x).translate()) or not len(x)%3==0 else str(Seq(x).translate()))
        df['productive'] = df['cdr3'].apply(lambda x: False if x=='None' else True)

        filter_barcode_rep = f'{self.outdir}/../03.assemble/{self.sample}_barcode_filter_report.tsv'
        filter_report = f'{self.outdir}/../03.assemble/{self.sample}_filter_report.tsv'
        df_for_clono, _ = Summarize.filter_cell(df, self.seqtype, filter_report, filter_barcode_rep)
        return df_for_clono

    @utils.add_log
    def annotation_process(self):
        df_for_clono = self.parse_contig()
        annotation_summary = get_vj_annot(df_for_clono, self.chains, self.paired_groups)
        for anno in annotation_summary:
            self.add_metric(anno['name'], anno['value'], anno.get('total'), anno.get('help_info'))

        df_clonotypes = self.parse_clonotype()
        title = 'Clonetypes'
        table_dict = self.get_table_dict(
            title = title,
            table_id = 'clonetypes',
            df_table = df_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']]
        )
        self.add_data(table_dict=table_dict)

        df_clonotypes['ClonotypeID'] = df_clonotypes['ClonotypeID'].astype("int")
        df_clonotypes.sort_values(by=['ClonotypeID'], inplace=True)
        Barplot = Bar_plot(df_bar=df_clonotypes).get_plotly_div()
        self.add_data(Barplot=Barplot)

    @utils.add_log
    def mapping_process(self):
        run_mapping(self.rds,self.contig,self.sample,self.outdir,self.assign_file)
        meta = pd.read_csv(glob.glob(f'{self.outdir}/{self.sample}_meta.csv')[0])
        metaTB = meta[meta['CellTypes'].isin(self.Celltype)]
        mappedmeta = meta[meta['Class']=='T/BCR']
        mappedmetaTB = mappedmeta[mappedmeta['CellTypes'].isin(self.Celltype)]
        
        # Transcriptome_cell_number = meta.shape[0]
        TB_cell_number = metaTB.shape[0]
        # Mapped_Transcriptome_cell_number = mappedmeta.shape[0]
        Mapped_TB_cell_number = mappedmetaTB.shape[0]
        mapping_summary=[]

        # mapping_summary.append({
        #     'item': 'Cell Number in Matched transcriptome',
        #     'count': Transcriptome_cell_number,
        #     'total_count': np.nan
        # })
        # mapping_summary.append({
        #     'item': 'Cell Number Successfully Mapped to transcriptome',
        #     'count': Mapped_Transcriptome_cell_number,
        #     'total_count': Transcriptome_cell_number
        # })
        
        mapping_summary.append({
            'item': f'{self._name} Number in Matched transcriptome',
            'count': TB_cell_number,
            'total_count': np.nan
        })
        mapping_summary.append({
            'item': f'Cell Number Successfully Mapped to {self._name} in transcriptome',
            'count': Mapped_TB_cell_number,
            'total_count': TB_cell_number
        })

        stat_file = self.outdir + '/mapping.txt'
        sum_df = pd.DataFrame(mapping_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file) 

    def run(self):
        self.annotation_process()
        
        try:
            if self.rds and self.assign_file:
                self.mapping_process()
        except AttributeError:
            print("rds file and type file do not exist" + "\n" )
        except ZeroDivisionError:
            print(f"Not found auto-assigned {self._name} in matched sc-RNA")


def annotation(args):
    with Annotation(args, display_title="V(D)J Annotation") as runner:
        runner.run()


def get_opts_annotation(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR',
                        choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
    return parser