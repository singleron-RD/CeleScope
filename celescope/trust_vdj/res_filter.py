import pandas as pd
from celescope.tools.Step import Step, s_common
from celescope.tools import utils


@utils.add_log
def beauty_res(outdir, barcode_report):
    res = pd.read_csv(barcode_report, sep='\t')
    rows = res.shape[0]
    loci = ['A', 'B']
    chians = ['chain2', 'chain1']
    for l in range(len(loci)):
        chain = chians[l]
        locus = loci[l]

        Vgenes, Dgenes, Jgenes, Cgenes, cdr3nts, cdr3aas, readcounts, fuls = [], [], [], [], [], [], [], []

        for i in range(rows):
            attr = res.loc[i, chain]
            attrs = attr.split(',')
            if len(attrs) == 10:
                V, D, J, C, cdr3nt, cdr3aa, readcount, fl = attrs[0], attrs[1], attrs[2], attrs[3], attrs[4], attrs[5], attrs[6], attrs[-1]
                Vgenes.append(V)
                Dgenes.append(D)
                Jgenes.append(J)
                Cgenes.append(C)
                cdr3nts.append(cdr3nt)
                cdr3aas.append(cdr3aa)
                readcounts.append(readcount)
                fuls.append(fl)
            elif len(attrs) != 10:
                Vgenes.append('NAN')
                Dgenes.append('NAN')
                Jgenes.append('NAN')
                Cgenes.append('NAN')
                cdr3nts.append('NAN')
                cdr3aas.append('NAN')
                readcounts.append('NAN')
                fuls.append('NAN')
            
        res[f'TR{locus}_V'] = Vgenes
        res[f'TR{locus}_D'] = Dgenes
        res[f'TR{locus}_J'] = Jgenes
        res[f'TR{locus}_C'] = Cgenes
        res[f'TR{locus}_cdr3nt'] = cdr3nts
        res[f'TR{locus}_cdr3aa'] = cdr3aas
        res[f'TR{locus}_readcount'] = readcounts
        res[f'TR{locus}_fl'] = fuls

    res.to_csv(f'{outdir}/new_barcode_report.tsv', sep='\t')

    return res


class Res_filter(Step):

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.sample = args.sample


    @utils.add_log
    def run(self):
        barcode_report = f'{self.outdir}/../02.trust_assemble/TRUST4/{self.sample}_barcode_report.tsv'
        res = beauty_res(self.outdir, barcode_report)
        filtered = res[(res['TRB_fl']!='0')&(res['TRA_fl']!='0')]
        fre = [''] * filtered.shape[0]
        filtered.insert(filtered.shape[1], 'Frequent', fre)

        clones = filtered.groupby(['TRA_cdr3aa', 'TRB_cdr3aa']).agg({'Frequent': 'count'})
        clones = clones.sort_values(by='Frequent', ascending=False)

        clones.to_csv(f'{self.outdir}/clonetype.tsv', sep='\t')


@utils.add_log
def res_filter(args):
    step_name = 'res_filter'
    res_filter_obj = Res_filter(args, step_name)
    res_filter_obj.run()


def get_opts_res_filter(parser, sub_program):
	if sub_program:
		parser = s_common(parser)