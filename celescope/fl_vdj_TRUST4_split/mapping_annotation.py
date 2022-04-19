import pandas as pd
import numpy as np 
import subprocess
import glob
from Bio.Seq import Seq
from celescope.fl_vdj_TRUST4_split.summarize import Summarize
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.fl_vdj_TRUST4_split.__init__ import TOOLS_DIR, CHAIN, PAIRED_CHAIN
from celescope.fl_vdj_TRUST4_split import trust_utils as tr
from celescope.tools.plotly_plot import Bar_plot


def run_mapping(rds, contig, sample, outdir, assign):
    cmd = (
        f'Rscript {TOOLS_DIR}/VDJmapping.R '
        f'--rds {rds} '
        f'--VDJ {contig} '
        f'--sample {sample} '
        f'--outdir {outdir} '
        f'--assign_file {assign}'
    )
    subprocess.check_call(cmd, shell=True)


class Mapping_annotation(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.match_dir = args.match_dir
        self.chains, self.paired_groups = Summarize._parse_seqtype(self.seqtype)
        self.min_read_count = args.min_read_count

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

    @staticmethod
    def _parse_seqtype(seqtype):
        return CHAIN[seqtype], PAIRED_CHAIN[seqtype]

    @staticmethod
    def parse_clonotype(outdir):
        df_clonotypes=pd.read_csv(f'{outdir}/../04.summarize/clonotypes.csv', sep=',')
        df_clonotypes['ClonotypeID'] = df_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        df_clonotypes['Frequency'] = df_clonotypes['frequency']
        df_clonotypes['Proportion'] = df_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        df_clonotypes['CDR3_aa'] = df_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))
        return df_clonotypes

    def parse_contig(self):
        df = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/{self.sample}_contig.csv', sep='\t', header=None) 
        df.columns = ['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis']
        df['d_gene'] = df['d_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['c_gene'] = df['c_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['cdr3'] = df['cdr3_nt'].apply(lambda x: 'None' if "*" in str(Seq(x).translate()) or not len(x)%3==0 else str(Seq(x).translate()))
        df['productive'] = df['cdr3'].apply(lambda x: False if x=='None' else True)

        filterbc_rep = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/barcoderepfl.tsv',sep='\t')
        filtered_report_out = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/trust_filter_report.out',sep='\t')
        df_for_clono = Summarize.filter_cell(df, self.seqtype, filterbc_rep, filtered_report_out)
        return df_for_clono

    @utils.add_log
    def annotation_process(self):
        df_for_clono = self.parse_contig()
        annotation_summary = tr.get_vj_annot(df_for_clono, self.chains, self.paired_groups)
        for anno in annotation_summary:
            self.add_metric(anno['name'], anno['value'], anno.get('total'), anno.get('help_info'))

        df_clonotypes = self.parse_clonotype(self.outdir)
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
        
        Transcriptome_cell_number = meta.shape[0]
        TB_cell_number = metaTB.shape[0]
        Mapped_Transcriptome_cell_number = mappedmeta.shape[0]
        Mapped_TB_cell_number = mappedmetaTB.shape[0]
        mapping_summary=[]

        mapping_summary.append({
            'item': 'Cell Number in Matched transcriptome',
            'count': Transcriptome_cell_number,
            'total_count': np.nan
        })
        mapping_summary.append({
            'item': 'Cell Number Successfully Mapped to transcriptome',
            'count': Mapped_Transcriptome_cell_number,
            'total_count': Transcriptome_cell_number
        })
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


def mapping_annotation(args):
    with Mapping_annotation(args, display_title="V(D)J Annotation") as runner:
        runner.run()

    
def get_opts_mapping_annotation(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR',
                        choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--min_read_count', help ='filter cell by read count number, int type required', default='auto')
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
    return parser