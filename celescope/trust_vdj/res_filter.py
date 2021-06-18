import pandas as pd
from celescope.tools.Step import Step, s_common
from celescope.tools import utils
from collections import defaultdict


@utils.add_log
def beauty_report(barcode_report):
    df = pd.read_csv(barcode_report, sep='\t')
    rows = df.shape[0]
    chains = ['chain2', 'chain1']
    dic = defaultdict(list)

    for l in range(len(chains)):
        chain = chains[l]

        items = {'V': 0, 'D': 1, 'J': 2, 'C': 3, 'CDR3nt': 4, 'CDR3aa': 5, 'readcount': 6, 'full_length_assembly': -1}        

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

    return res


def get_clone_table(df, Seqtype):
    res = pd.DataFrame()
    group_type = []
    if Seqtype == 'TCR':
        chains = ['TRA', 'TRB']
    if Seqtype == 'BCR':
        chains = ['IGH', 'IGL', 'IGK']
        for chain in chains:
            tmp = df[df['V'].str.contains(chain, na=False)]
            tmp = tmp.set_index('barcode')
            tmp = tmp.rename(columns=lambda x: f'{chain}_'+x)

            res = pd.concat([res, tmp], axis=1, join='outer', sort=False).fillna('None')
            group_type.append(f'{chain}_CDR3aa')

    
    Frequent = [''] * res.shape[0]
    res.insert(res.shape[1], 'Frequent', Frequent)
    clonetypes = res.groupby(group_type).agg({'Frequent': 'count'})
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
    
    clonetypes['CloneId'] = [i for i in range(1, (clonetypes.shape[0]+1))]
    clonetypes = clonetypes.reindex(columns=list(['CloneId', 'TRA_CDR3aa', 'TRB_CDR3aa', 'Frequent', 'Proportion']))   

    return clonetypes


class Res_filter(Step):

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.sample = args.sample
        self.Seqtype = args.Seqtype


    @utils.add_log
    def run(self):
        barcode_report = f'{self.outdir}/../02.trust_assemble/TRUST4/{self.sample}_barcode_report.tsv'
        df = beauty_report(barcode_report)
        df.to_csv(f'{self.outdir}/{self.sample}_barcode_report.tsv', sep='\t')

        clones = get_clone_table(df, self.Seqtype)

        clones.to_csv(f'{self.outdir}/clonetype.tsv', sep='\t')

        title = 'Clonetypes'
        table_dict = self.get_table(title, 'clonetypes_table', clones)

        self.add_data_item(table_dict=table_dict)

        self.clean_up()


@utils.add_log
def res_filter(args):
    step_name = 'res_filter'
    res_filter_obj = Res_filter(args, step_name)
    res_filter_obj.run()


def get_opts_res_filter(parser, sub_program):
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        parser = s_common(parser)