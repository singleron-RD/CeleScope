import pandas as pd
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.tools.cellranger3 import get_plot_elements
import numpy as np
import pysam
import subprocess
import os


def trans_rep(report, outprefix):
    cmd = (
        f'perl /SGRNJ03/randd/zhouxin/software/TRUST4/trust-barcoderep-to-10X.pl '
        f'{report} '
        f'{outprefix}'
    )
    subprocess.check_call(cmd, shell=True)


class Res_format(Step):
    """
    Features

    - Calculate clonetypes.

    Output
    - `05.res_format/clonetypes.tsv` Record each clonetype and its frequent.
    - `05.res_format/{sample}_barcode_report.tsv` Record detailed chain information of each barcode.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.sample = args.sample
        self.Seqtype = args.Seqtype
        self.full_length = args.full_length
        self.report = args.report
        self.fa = args.fa
        self.count_file = args.count_file
        self.cdr3out = args.cdr3out

        if self.Seqtype == 'TCR':
            self.string = 't'
            self.chain = ['TRA', 'TRB']
            self.paired_groups = ['TRA_TRB']
        elif self.Seqtype == 'BCR':
            self.string = 'b'
            self.chain = ['IGH', 'IGL', 'IGK']
            self.paired_groups = ['IGH_IGL', 'IGH_IGK']

        self.all_out_rep_prefix = f'{self.outdir}/all'
        self.filter_out_rep_prefix  = f'{self.outdir}/filter'
        self.all_rep = f'{self.all_out_rep_prefix}_{self.string}.csv'
        self.filter_rep = f'{self.filter_out_rep_prefix}_{self.string}.csv'

    @utils.add_log
    def get_all_rep(self):
        trans_rep(self.report, self.all_out_rep_prefix)


    @utils.add_log
    def get_filter_rep(self):
        out_rep = f'{self.outdir}/{self.sample}_filter_report.tsv'
        cmd = (
            f'perl /SGRNJ03/randd/zhouxin/software/TRUST4/trust-barcoderep.pl '
            f'{self.cdr3out} '
            f'-a {self.fa} '
            f'--noImputation > '
            f'{out_rep}'
        )
        subprocess.check_call(cmd, shell=True)
        trans_rep(out_rep, self.filter_out_rep_prefix)
        os.remove(out_rep)


    @utils.add_log
    def get_len(self):
        with pysam.FastaFile(self.fa) as fh:
            res = {}
            names = fh.references
            lengths = fh.lengths
            res['contig_id'] = names
            res['length'] = lengths
            
            df = pd.DataFrame(res, columns=list(res.keys()))
            
            data = pd.read_csv(self.filter_rep, sep=',')
            data = data.drop('length', 1)
            df_len = pd.merge(df, data, on='contig_id', how='inner')
            df_len.to_csv(self.filter_rep, sep=',', index=False)
            
            all_res = pd.read_csv(self.all_rep, sep=',')
            all_res = all_res.drop('length', 1)
            df_all_len = pd.merge(df, all_res, on='contig_id', how='inner')
            df_all_len.to_csv(self.all_rep, sep=',', index=False)
            

    @utils.add_log
    def get_clonetypes(self):
        data = pd.read_csv(self.filter_rep, sep=',')
        data = data[['barcode', 'chain', 'cdr3', 'cdr3_nt']]
        with open(f'{self.outdir}/clonetypes.csv', 'w') as fh:
            fh.write('barcode,cdr3s_aa,cdr3s_nt\n')
            barcodes = set(data['barcode'].tolist())
            data = data.set_index('barcode', drop=False)
            for cb in barcodes:
                aa = []
                nt = []
                for c in self.chain:
                    tmp = data[(data['barcode']==cb)&(data['chain']==c)]
                    if not tmp.empty:
                        cdr3_aa = tmp.loc[cb, 'cdr3']
                        cdr3_nt = tmp.loc[cb, 'cdr3_nt']
                        aa.append(f'{c}:{cdr3_aa}')
                        nt.append(f'{c}:{cdr3_nt}')
                aas = ';'.join(aa)
                nts = ';'.join(nt)
                string = f'{cb},{aas},{nts}'
                fh.write(string+'\n')
                
            fh.close()

        df = pd.read_csv(f'{self.outdir}/clonetypes.csv', sep=',')
        df = df.groupby(['cdr3s_aa', 'cdr3s_nt'], as_index=False).agg({'barcode': 'count'})
        df = df.rename(columns={'barcode': 'frequency'})
        df = df.sort_values(by='frequency', ascending=False)
        sum_c = df['frequency'].sum()
        df['proportion'] = df['frequency'].apply(lambda x: x/sum_c)
        df['clonotype_id'] = [f'clonetype{i}' for i in range(1, df.shape[0]+1)]

        df = df.reindex(columns=['clonotype_id', 'frequency', 'proportion', 'cdr3s_aa', 'cdr3s_nt'])
        df.to_csv(f'{self.outdir}/clonetypes.csv', sep=',', index=False)


    @utils.add_log
    def gen_table(self):
        df = pd.read_csv(f'{self.outdir}/clonetypes.csv', sep=',')
        df_table = pd.DataFrame()
        def function_x(x):
            return x.split(';')
        df_table['Clonotype ID'] = df['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        df_table['Chain1'] = df['cdr3s_aa'].apply(lambda x: f'{function_x(x)[0]}')
        df_table['Chain2'] = df['cdr3s_aa'].apply(lambda x: f'{function_x(x)[1]}' if (len(function_x(x))==2) else 'None')
        df_table['Frequency'] = df['frequency']
        df_table['Proportion'] = df['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        title = 'Clonetypes'

        table_dict = self.get_table(title, 'clonetypes_table', df_table)

        self.add_data_item(table_dict=table_dict)


    @utils.add_log
    def gen_stat_file(self):
        res_format_summary = []
        data = pd.read_csv(self.filter_rep, sep=',')
        data = data[['barcode', 'chain']]
        barcodes = set(data['barcode'].tolist())
        total_count = len(barcodes)

        res = pd.DataFrame()
        for c in self.chain:
            tmp = data[data['chain']==c]
            item = f'Number of cells with {c}'
            count = int(tmp.shape[0])
            res_format_summary.append({'item': item, 'count': count, 'total_count': total_count})

            tmp = tmp.set_index('barcode')
            tmp = tmp.rename(columns={'chain': f'{c}_chain'})

            res = pd.concat([res, tmp], axis=1, join='outer', sort=False).fillna('None')

        for pg in self.paired_groups:
            attrs = pg.split('_')
            chain1 = attrs[0]
            chain2 = attrs[1]
            tmp = res[(res[f'{chain1}_chain']!='None') & (res[f'{chain2}_chain']!='None')]
            item = f'Number of cells with paired {chain1} and {chain2}'
            count = int(tmp.shape[0])
            res_format_summary.append({
                'item': item,
                'count': count,
                'total_count': total_count
            })

        stat_file = self.outdir + '/stat.txt'

        sum_df = pd.DataFrame(res_format_summary, columns=['item', 'count', 'total_count'])

        utils.gen_stat(sum_df, stat_file)

        #plot

        df_umi = pd.read_csv(self.count_file, sep='\t')
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in barcodes else 'UB')
        df_umi.to_csv(self.count_file, sep='\t', index=False)
        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(self.count_file))


    @utils.add_log
    def run(self):
        self.get_all_rep()
        self.get_filter_rep()
        self.get_len()
        self.get_clonetypes()
        self.gen_table()
        self.gen_stat_file()

        self.clean_up()


@utils.add_log
def res_format(args):
    step_name = 'res_format'
    res_format_obj = Res_format(args, step_name)
    res_format_obj.run()


def get_opts_res_format(parser, sub_program):
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--full_length', help='only output full length assembly', action='store_true')
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--report', help='assemble report', required=True)
        parser.add_argument('--fa', help='assembled fasta file', required=True)
        parser.add_argument('--count_file', help='UMI count file', required=True)
        parser.add_argument('--cdr3out', help='cdr3 out file', required=True)






    
