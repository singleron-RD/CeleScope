
import pandas as pd
from celescope.tools import utils
from celescope.tools.step import Step, s_common

DIR = '/SGRNJ03/randd/zhouxin/software/TRUST4/'


class Res_sum(Step):
    """
    ## Features

    - Calculate clonetypes.

    ## Output
    - `06.res_sum/clonetypes.tsv` Record each clonetype and its frequent.
    - `06.res_sum/all_{type}.csv` Containing detailed information for each barcode.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.sample = args.sample
        self.Seqtype = args.Seqtype
        self.all_rep = args.all_rep
        self.fa = args.fa

        if self.Seqtype == 'TCR':
            self.string = 't'
            self.chain = ['TRA', 'TRB']
            self.paired_groups = ['TRA_TRB']
        elif self.Seqtype == 'BCR':
            self.string = 'b'
            self.chain = ['IGH', 'IGL', 'IGK']
            self.paired_groups = ['IGH_IGL', 'IGH_IGK']

    @utils.add_log
    def gen_stat_file(self):
        res_sum_summary = []
        data = pd.read_csv(self.all_rep, sep=',')
        pro_data = data[(data['full_length'] == True) & (data['productive'] == True)]
        barcodes = set(data['barcode'].tolist())
        total_count = len(barcodes)
        res = pd.DataFrame()

        for c in self.chain:
            tmp = pro_data[pro_data['chain'] == c]

            tmp = tmp.set_index('barcode')
            tmp = tmp.rename(columns={'chain': f'{c}_chain'})

            res = pd.concat([res, tmp], axis=1, join='outer', sort=False).fillna('None')

        for pg in self.paired_groups:
            attrs = pg.split('_')
            chain1 = attrs[0]
            chain2 = attrs[1]
            tmp = res[(res[f'{chain1}_chain'] != 'None') & (res[f'{chain2}_chain'] != 'None')]
            item = f'Number of cells with V-J spanning productive paired {chain1} and {chain2}'
            count = int(tmp.shape[0])
            res_sum_summary.append({
                'item': item,
                'count': count,
                'total_count': total_count
            })

        for c in self.chain:
            dic = {}
            tmp = data[data['chain'] == c]
            dic[f'Cells With {c} Contig'] = tmp.shape[0]
            dic[f'Cells With V-J spanning {c} Contig'] = tmp[tmp['full_length'] == True].shape[0]
            dic[f'Cells With productive {c} Contig'] = tmp[(
                tmp['productive'] == True) & (tmp['full_length'] == True)].shape[0]

            for item in list(dic.keys()):
                res_sum_summary.append({
                    'item': item,
                    'count': int(dic[item]),
                    'total_count': total_count
                })

        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(res_sum_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)

    @utils.add_log
    def run(self):

        self.gen_stat_file()

        self._clean_up()


@utils.add_log
def res_sum(args):
    res_sum_obj = Res_sum(args)
    res_sum_obj.run()


def get_opts_res_sum(parser, sub_program):
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--all_rep', help='filtered assemble report without imputation', required=True)
        parser.add_argument('--fa', help='assembled fasta file', required=True)
