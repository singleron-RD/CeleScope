import pandas as pd
from celescope.tools import utils
from celescope.tools import count as super_count
from celescope.tools.matrix import ROW

class Count(super_count.Count):
    """
    ## Features
    - Generate expression matrix.
    - Well statstic.

    ## Output
    - `{sample}_raw_feature_bc_matrix` The expression matrix of all detected barcodes in [Matrix Market Exchange Formats](
        https://math.nist.gov/MatrixMarket/formats.html). 
    - `{sample}_count_detail.txt.gz` 4 columns: 
        - barcode  
        - gene ID  
        - UMI count  
        - read_count  
    - `{sample}_counts.txt` 6 columns:
        - Barcode: barcode sequence
        - read: read count of each barcode
        - UMI: UMI count for each barcode
        - geneID: gene count for each barcode
        - mark: cell barcode or backgound barcode.
            `CB` cell  
            `UB` background  
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        self.raw_count_file = f'{self.outdir}/{self.sample}_counts.txt'
        self.marked_count_file = f'{self.outdir}/{self.sample}_counts_report.txt'
        self.umi_cutoff = args.umi_cutoff
        self.read_cutoff = args.read_cutoff
        self.gene_cutoff = args.gene_cutoff
        self._table_id = self.assay

    
    @utils.add_log
    def run(self):
        ## output exprssion matrix
        df = pd.read_table(self.args.count_detail, header=0, dtype={'geneID':str}, index_col=[0,1])
        self.write_sparse_matrix(df, self.raw_matrix_dir)
        ## output stats
        df_bc = self.get_df_bc(df)
        df_bc.index.name = 'Well'
        sort_col = ['UMI', 'read', ROW]
        df_bc = df_bc[sort_col]
        df_bc.columns = ['UMI', 'read', 'gene']
        df_bc.to_csv(self.raw_count_file, sep='\t')
        df_bc_valid = df_bc[ (df_bc['UMI']>=self.umi_cutoff) & (df_bc['read']>=self.read_cutoff) & (df_bc['gene']>=self.gene_cutoff) ]
        df_bc_valid.to_csv(self.marked_count_file, sep='\t')
        ## report
        self.add_count_metrics(df_bc_valid)

    @utils.add_log
    def add_count_metrics(self, df):
        out_condi = ' '
        out_condi += f'umi>={self.umi_cutoff} ' if self.umi_cutoff>0 else ''
        out_condi += f'read>={self.read_cutoff} ' if self.read_cutoff>0 else ''
        out_condi += f'gene>={self.gene_cutoff} ' if self.gene_cutoff>0 else ''
        self.add_help_content(
            name='Well:', content=f'The 1st column in the table is the well barcode sequence. Only output wells that meet the conditions({out_condi}).')
        stats = df.describe()
        stats.columns = ['UMI', 'Reads', 'Genes']
        for item in ['Reads', 'UMI', 'Genes']:
            self.add_metric(
                name=f'Median {item} per Well',
                value=int(stats.loc['50%',item]),
                help_info=''
            )
        for item in ['Reads', 'UMI', 'Genes']:
            self.add_metric(
                name=f'Mean {item} per Well',
                value=int(stats.loc['mean',item]),
                help_info=''
            )
        table_dict = self.get_table_dict(
            title='Detailed information per Well',
            table_id=self._table_id,
            df_table=df.reset_index(),
        )
        self.add_data(table_dict=table_dict)



def count(args):
    with Count(args, display_title="Wells") as runner:
        runner.run()

def get_opts_count(parser, sub_program):
    parser.add_argument('--umi_cutoff', default=500, type=int,
        help='If the UMI number exceeds the threshold, it is considered a valid well and reported.'
    )
    parser.add_argument('--gene_cutoff', default=0, type=int,
        help='If the gene number exceeds the threshold, it is considered a valid well and reported.'
    )
    parser.add_argument('--read_cutoff', default=0, type=int,
        help='If the read number exceeds the threshold, it is considered a valid well and reported.'
    )
    super_count.get_opts_count(parser, sub_program)

