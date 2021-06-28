import pandas as pd
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from collections import defaultdict
from celescope.tools.cellranger3 import get_plot_elements
import numpy as np
import pysam


def get_len(fa):
    with pysam.FastaFile(fa) as fh:
        res = {}
        names = fh.references
        lengths = fh.lengths
        res['contig_id'] = names
        res['length'] = lengths
        
        df = pd.DataFrame(res, columns=list(res.keys()))
        return df


@utils.add_log
def beauty_report(barcode_report, fa):
    df_len = get_len(fa)
    df = pd.read_csv(barcode_report, sep='\t')
    rows = df.shape[0]
    chains = ['chain2', 'chain1']
    dic = defaultdict(list)

    for l in range(len(chains)):
        chain = chains[l]

        items = {'V': 0, 'D': 1, 'J': 2, 'C': 3, 'CDR3nt': 4, 'CDR3aa': 5, 'readcount': 6, 'contig_id': -3, 'full_length_assembly': -1}        

        for i in range(rows):
            cb = df.loc[i, '#barcode']
            dic['barcode'].append(cb)
            for item in items:
                attr = df.loc[i, chain]
                attrs = attr.split(',')

                if len(attrs) == 10:
                    dic[f'{item}'].append(attrs[items[item]])

                elif len(attrs) != 10:
                    dic[f'{item}'].append('None')

    res = pd.DataFrame(dic, columns=list(dic.keys()))

    df_res = pd.merge(res, df_len, on='contig_id', how='inner')

    return df_res

@utils.add_log
def get_clone_table(df, Seqtype):
    res_filter_summary = []

    res = pd.DataFrame()
    group_type = []
    if Seqtype == 'TCR':
        chains = ['TRA', 'TRB']
        paired_groups = ['TRA_TRB']
    if Seqtype == 'BCR':
        chains = ['IGH', 'IGL', 'IGK']
        paired_groups = ['IGH_IGL', 'IGH_IGK']
    for c in chains:
        tmp = df[df['V'].str.contains(c, na=False)]
        tmp = tmp.set_index('barcode')
        tmp = tmp.rename(columns=lambda x: f'{c}_{x}')

        res = pd.concat([res, tmp], axis=1, join='outer', sort=False).fillna('None')
        group_type.append(f'{c}_CDR3aa')
    
    Frequent = [''] * res.shape[0]
    res.insert(res.shape[1], 'Frequent', Frequent)
    clonetypes = res.groupby(group_type, as_index=False).agg({'Frequent': 'count'})
    clonetypes = clonetypes.sort_values(by='Frequent', ascending=False)

    sum_c = clonetypes['Frequent'].sum()
    proportions = []
    for f in list(clonetypes['Frequent']):
        p = f/sum_c
        p = p * 100
        p = round(p, 2)
        p = str(p) + '%'
        proportions.append(p)
    clonetypes['Proportion'] = proportions
    clonetypes = clonetypes.sort_values(by='Frequent', ascending=False)
    clonetypes = clonetypes.reset_index()
    
    clonetype_ids = [(i+1) for i in clonetypes.index.tolist()]
    clonetypes['index'] = clonetype_ids
    clonetypes = clonetypes.rename(columns={'index': 'CloneId'})  

    total_count = int(clonetypes['Frequent'].sum())

    res_filter_summary.append({
        'item': 'Estimated Number of Cells',
        'count': total_count,
        'total_count': np.nan
    })
    
    for group in group_type:
        chain = group.strip('_CDR3aa')
        tmp = clonetypes[clonetypes[group]!='None']
        count = int(tmp['Frequent'].sum())
        item = f'Cells with {chain}'
        res_filter_summary.append({
            'item': item,
            'count': count,
            'total_count': total_count
        })

    for pg in paired_groups:
        attrs = pg.split('_')
        chain1 = attrs[0]
        chain2 = attrs[1]
        tmp = clonetypes[(clonetypes[f'{chain1}_CDR3aa']!='None') & (clonetypes[f'{chain2}_CDR3aa']!='None')]
        item = f'Cells with paired {chain1} and {chain2}'
        count = int(tmp['Frequent'].sum())
        res_filter_summary.append({
            'item': item,
            'count': count,
            'total_count': total_count
        })


    return clonetypes, res_filter_summary


class Res_filter(Step):
    """
    Features

    - Calculate clonetypes.

    Output
    - `05.res_filter/clonetypes.tsv` Record each clonetype and its frequent.
    - `05.res_filter/{sample}_barcode_report.tsv` Record detailed chain information of each barcode.
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


    @utils.add_log
    def run(self):
        barcode_report = self.report
        fa = self.fa
        df = beauty_report(barcode_report, fa)

        if self.full_length:
            df = df[df['full_length_assembly']=='1']
        df.to_csv(f'{self.outdir}/{self.sample}_barcode_report.tsv', sep='\t', index=False)

        clones, res_filter_summary = get_clone_table(df, self.Seqtype)

        # plot barcode umi
        count_file = self.count_file
        df_umi = pd.read_csv(count_file, sep='\t', index_col=False)
        cells = set(df['barcode'].tolist())
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if (x in cells) else 'UB')
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi.to_csv(count_file, sep='\t', index=False)

        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_file))        

        if self.Seqtype == 'TCR':
            chains = ['TRA', 'TRB']
        elif self.Seqtype == 'BCR':
            chains = ['IGH', 'IGL', 'IGK']

        for chain in chains:
            tmp = df[df['V'].str.contains(chain, na=False)]
            barcodes = tmp['barcode'].tolist()
            if len(barcodes) != 0:
                df_bc = pd.DataFrame(barcodes, columns=['barcode'])
            else:
                continue

            tmp_df = pd.merge(df_umi, df_bc, on='barcode', how='inner')

            mid = int(tmp_df['UMI'].median())
            item = f'Median {chain} UMIs per cell'
            res_filter_summary.append({
                'item': item,
                'count': mid,
                'total_count': np.nan
            })

        match_file = open(f'{self.outdir}/{self.sample}_matched_barcodes.txt', 'r')
        matched_cells = match_file.readlines()
        res_filter_summary.insert(0, {'item': 'Number of matched cells', 
                                    'count': len(matched_cells),
                                    'total_count': np.nan
                                })

        clones.to_csv(f'{self.outdir}/clonetype.tsv', sep='\t')

        title = 'Clonetypes'
        table_dict = self.get_table(title, 'clonetypes_table', clones)

        self.add_data_item(table_dict=table_dict)


        stat_file = self.outdir + '/stat.txt'

        sum_df = pd.DataFrame(res_filter_summary, columns=['item', 'count', 'total_count'])

        utils.gen_stat(sum_df, stat_file)

        self.clean_up()


@utils.add_log
def res_filter(args):
    step_name = 'res_filter'
    res_filter_obj = Res_filter(args, step_name)
    res_filter_obj.run()


def get_opts_res_filter(parser, sub_program):
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--full_length', help='only output full length assembly', action='store_true')
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--report', help='assemble report', required=True)
        parser.add_argument('--fa', help='assembled fasta file', required=True)
        parser.add_argument('--count_file', help='UMI count file', required=True)






    
