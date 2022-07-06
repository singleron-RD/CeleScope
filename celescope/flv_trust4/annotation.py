import pandas as pd
import subprocess
import glob

from celescope.flv_trust4.summarize import Summarize
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.flv_trust4.__init__ import TOOLS_DIR
from celescope.tools.plotly_plot import Bar_plot


CELL_TYPE_DICT = {
    'TCR':['T_cells', 'NKT_cells', 'T cells', 'NK T cells', 'Tcells'],
    'BCR':['Plasma_cells', 'B_cells', 'Mature_B_cell', 'Plasma cells', 'B cells', 'Bcells'],
    }


class Annotation(Step):
    """
    ## Features

    - Output assembled T/B cells mapping to transcriptome if rds and auto-assign info exist in match directory.
    - Generate VDJ annotation info, clonetypes table and bar-plot of clonetypes distribution in html.

    ## Output
    - `05.annotation/{sample}_assign.png` Umap plot of Auto-assigned celltype in transcriptome.
    - `05.annotation/{sample}_cluster_umap.png` Umap plot of Cluster in transcriptome.
    - `05.annotation/{sample}_umapplot.png` Umap plot of assembled barcodes marked as read color.
    - `05.annotation/{sample}_distribution.txt` Number of assembled barcodes in every clusters.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.match_dir = args.match_dir
        self.chains, self.paired_groups = Summarize._parse_seqtype(self.seqtype)
        self.contig_file = f'{args.summarize_out}/{self.sample}_filtered_contig.csv'
        self.clonotype_file = f'{args.summarize_out}/clonotypes.csv'        
        self.celltype_set = CELL_TYPE_DICT[self.seqtype]

        try:
            self.rds = glob.glob(f'{self.match_dir}/06.analysis/*.rds')[0]
            self.assign_file = glob.glob(f'{self.match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')[0]
        except IndexError:
            pass

    def parse_clonotype(self):
        """Generate clonotypes table in html.
        """
        df_clonotypes=pd.read_csv(self.clonotype_file, sep=',')
        df_clonotypes['ClonotypeID'] = df_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        df_clonotypes['Frequency'] = df_clonotypes['frequency']
        df_clonotypes['Proportion'] = df_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        df_clonotypes['CDR3_aa'] = df_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))

        return df_clonotypes

    @staticmethod
    @utils.add_log
    def get_vdj_metric(df, chains, pairs):
        """
        Add vdj metrics in html.
        """
        metric_result = []
        fl_pro_pair_df = pd.DataFrame(df[df['productive']==True].barcode.value_counts())
        fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']>=2]
        cell_nums = len(set(df['barcode']))

        metric_result.append({
            'name': 'Cells With Productive V-J Spanning Pair',
            'value': fl_pro_pair_df.shape[0],
            'total': cell_nums,
        })

        for pair in pairs:
            chain1, chain2 = pair.split('_')[0], pair.split('_')[1]
            cbs1 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain1)].barcode)
            cbs2 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain2)].barcode)
            paired_cbs = len(cbs1.intersection(cbs2))

            metric_result.append({
                'name': f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
                'value': paired_cbs,
                'total': cell_nums,
                'help_info': "Fraction of cell-associated barcodes with one productive contig for each chain of the receptor pair.A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
            })

        for chain in chains:
        
            metric_result.append({
                'name': f'Cells With {chain} Contig',
                'value': len(set(df[df['chain']==chain].barcode)),
                'total': cell_nums,
                'help_info': f'Fraction of cell-associated barcodes with at least one {chain} contig annotated as a full or partial V(D)J gene'
            })
            metric_result.append({
                'name': f'Cells With CDR3-annotated {chain} Contig',
                'value': len(set(df[(df['chain']==chain)&(df['cdr3']!=None)].barcode)),
                'total': cell_nums,
            })
            metric_result.append({
                'name': f'Cells With V-J Spanning {chain} Contig',
                'value': len(set(df[(df['full_length']==True)&(df['chain']==chain)].barcode)),
                'total': cell_nums,
                'help_info': f"Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for {chain}"
            })
            metric_result.append({
                'name': f'Cells With Productive {chain} Contig',
                'value': len(set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain)].barcode)),
                'total': cell_nums,
                'help_info': "Fraction of cell-associated barcodes with productive IGL chain. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
            })

        return metric_result 

    @utils.add_log
    def annotation_process(self):
        """Add metrics, clonotypes table, and bar-plot of clonotypes distribution in html.
        """
        df = pd.read_csv(self.contig_file)

        annotation_summary = Annotation.get_vdj_metric(df, self.chains, self.paired_groups)
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

    @staticmethod
    @utils.add_log
    def run_mapping(rds, contig, sample, outdir, assign):
        cmd = (
            f'Rscript {TOOLS_DIR}/VDJmapping.R '
            f'--rds {rds} '
            f'--VDJ {contig} '
            f'--sample {sample} '
            f'--outdir {outdir} '
            f'--assign_file {assign} '
            '2>&1 '
        )
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def mapping_process(self):
        """Mapping result with matched scRNA.
        """
        Annotation.run_mapping(
            self.rds, self.contig_file, self.sample, self.outdir, self.assign_file
            )

        meta = pd.read_csv(f'{self.outdir}/{self.sample}_meta.csv')
        metaTB = meta[meta['CellTypes'].isin(self.celltype_set)]
        mappedmeta = meta[meta['Class']=='T/BCR']
        mappedmetaTB = mappedmeta[mappedmeta['CellTypes'].isin(self.celltype_set)]
        
        self.add_metric(
            'Total Cell Number in Matched transcriptome',
            meta.shape[0],
            show=False)
        self.add_metric(
            'Cell Number Successfully Mapped to transcriptome',
            mappedmeta.shape[0],
            show=False)
        self.add_metric(
            'T/B cell Number in Matched transcriptome',
            metaTB.shape[0],
            show=False)
        self.add_metric(
            'Cell Number Successfully Mapped to T/B cell in transcriptome',
            mappedmetaTB.shape[0],
            show=False)
    
    def run(self):
        self.annotation_process()

        try:
            if self.rds and self.assign_file:
                self.mapping_process()
        except AttributeError:
            print("rds file and type file do not exist" + "\n" )
        except ZeroDivisionError:
            print("Not found auto-assigned T/B cell in matched sc-RNA")


@utils.add_log
def annotation(args):
    with Annotation(args, display_title="V(D)J Annotation") as runner:
        runner.run()


def get_opts_annotation(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR',
                        choices=['TCR', 'BCR'], required=True)               
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
        parser.add_argument('--summarize_out', help='summarize output directory', required=True)  
    return parser