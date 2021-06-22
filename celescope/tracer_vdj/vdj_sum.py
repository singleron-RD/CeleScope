import os
from re import I
import pandas as pd
from Bio.Seq import Seq
import numpy as np
from celescope.tools import utils
from celescope.tools.Step import Step, s_common
from celescope.tools.cellranger3 import get_plot_elements
import glob
import pysam



def get_umi_count(fq):
    umis = []
    with pysam.FastxFile(fq) as fh:
        for entry in fh:
            name = entry.name
            name = name.split('_')
            umi = name[1]
            umis.append(umi)
    count = len(set(umis))
             
    return count


@utils.add_log
def tpm_count(ass_dir):
    rec = pd.read_csv(f'{ass_dir}/tracer/filtered_TCRAB_summary/recombinants.txt', sep='\t')  
    # ass_dir outdir/sample/04.go_assemble
    productive = rec[rec['productive'] == True]
    indx = list(productive.index)
    tpms = []
    for i in indx:
        cell_name = productive.loc[i, 'cell_name']
        rec_id = productive.loc[i, 'recombinant_id']
        with open(f'{ass_dir}/tracer/{cell_name}/expression_quantification/abundance.tsv') as tsvf:
            for line in tsvf:
                if rec_id in line:
                    line = line.rstrip()
                    line = line.split('\t')
                    tpm = float(line[4])
                    tpms.append(tpm)
    productive.insert(loc=productive.shape[1], column='TPM', value=tpms)

    return productive


@utils.add_log
def filtering(Seqtype, ass_dir, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if Seqtype == 'TCR':
        data = tpm_count(ass_dir)
        cell_name = sorted(list(set(list(data['cell_name']))))
        filtered = pd.DataFrame()
        df = pd.DataFrame(cell_name, columns=['cell_name'])
        loci = ['A', 'B']
        for locus in loci:
            tmp = data[data['locus']==locus]
            tmp = tmp.sort_values(by='TPM', ascending=False)
            tmp = tmp.drop_duplicates('cell_name', 'first')
            filtered = filtered.append(tmp, ignore_index=True)

            tmp = tmp.rename(columns={'CDR3aa': f'TR{locus}_CDR3aa'})
            clones = tmp[['cell_name', f'TR{locus}_CDR3aa']]
            df = pd.merge(df, clones, on='cell_name', how='outer')

        df = df.fillna('None')

        clonetypes = df.groupby(['TRA_CDR3aa', 'TRB_CDR3aa']).agg({'cell_name': 'count'})
        clonetypes = clonetypes.sort_values(by='cell_name', ascending=False)
        clonetypes = clonetypes.rename(columns={'cell_name': 'Frequency'})

        clonetypes.to_csv(f'{outdir}/clonetypes.tsv', sep='\t')
        filtered.to_csv(f'{outdir}/filtered.txt', sep='\t')

    elif Seqtype == 'BCR':

        data = pd.read_csv(f'{ass_dir}/bracer/filtered_BCR_summary/changeodb.tab', sep='\t')
        data = data[(data['FUNCTIONAL'] == True) & (data['IN_FRAME'] == True)]
        cell_name = sorted(list(set(data['CELL'].tolist())))
        filtered = pd.DataFrame()

        tmp = data[data['LOCUS'] == 'H']
        tmp = tmp.sort_values(by='TPM', ascending=False)
        tmp = tmp.drop_duplicates('CELL', 'first')
        filtered = filtered.append(tmp, ignore_index=True)

        tmp2 = data[data['LOCUS'] != 'H']
        tmp2 = tmp2.sort_values(by='TPM', ascending=False)
        tmp2 = tmp2.drop_duplicates('CELL', 'first')
        filtered = filtered.append(tmp2, ignore_index=True)

        df = pd.DataFrame(cell_name, columns=['CELL'])

        loci = ['H', 'L', 'K']
        for locus in loci:
            tmp = filtered[filtered['LOCUS'] == locus][['CELL', 'JUNCTION']]
            tmp.columns = ['CELL', f'JUNCTION_{locus}']
            ntseqs = tmp[f'JUNCTION_{locus}'].tolist()
            tmplist = []
            for nt in ntseqs:
                nt = Seq(nt)
                nt = nt.translate()
                tmplist.append(str(nt))
            tmp.insert(tmp.shape[1], f'IG{locus}_CDR3aa', tmplist)

            df = pd.merge(df, tmp, on='CELL', how='outer')

        df = df.fillna('None')
            
        clonetypes = df.groupby(['IGH_CDR3aa', 'IGL_CDR3aa', 'IGK_CDR3aa']).agg({'CELL': 'count'})
        clonetypes = clonetypes.sort_values(by='CELL', ascending=False)
        clonetypes = clonetypes.rename(columns={'CELL': 'Frequency'})

        clonetypes.to_csv(f'{outdir}/clonetypes.tsv', sep='\t')
        filtered = filtered.rename(columns={'CELL': 'cell_name'})
        filtered.to_csv(f'{outdir}/filtered.txt', sep='\t')


class Vdj_sum(Step):
    """
    Features

    - Filter tracer results by TPM.
    - Calculate clonetypes.

    Output

    - `05.vdj_sum/filtered.txt` Filtered results of tracer. Each cell has unique chain for each locus.
    - `05.vdj_sum/clonetypes.txt` Clonetypes calculation. 5 (TCR) or 6 (BCR) columns, clonetypeId, (detailed clonetypes), frequency, proportion.
    """
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        self.Seqtype = args.Seqtype
        self.fastq_dir = args.fastq_dir
        self.ass_dir = args.ass_dir


    @utils.add_log                    
    def run(self):
        ass_dir = self.ass_dir
        outdir = self.outdir
        fastq_dir = self.fastq_dir
        Seqtype = self.Seqtype

        filtering(Seqtype, ass_dir, outdir)

        filter_data = pd.read_csv(f'{outdir}/filtered.txt', sep='\t')

        stat_file = outdir + '/stat.txt'

        vdj_sum_summary = []
        
        count_umi_file = f'{fastq_dir}/../{self.sample}_count.txt'
        count_umi = pd.read_csv(count_umi_file, sep='\t', index_col=0)
        median_all = int(count_umi['UMI'].median())

        clonetypes = pd.read_csv(f'{outdir}/clonetypes.tsv', sep='\t')

        productive_cells = set(filter_data['cell_name'].tolist())
        productive_cells_num = len(productive_cells)

        if Seqtype == 'TCR':
            # barcode umi plot
            count_umi['mark'] = count_umi['cell_name'].apply(lambda x: "CB" if (x in productive_cells) else "UB")

            count_umi.to_csv(count_umi_file, sep='\t')

            self.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_umi_file))

            # clonetype table
            sum_c = clonetypes['Frequency'].sum()
            proportions = []
            for f in list(clonetypes['Frequency']):
                p = f/sum_c
                p = p * 100
                p = round(p, 2)
                p = str(p) + '%'
                proportions.append(p)
            clonetypes['Proportion'] = proportions
            clonetypes = clonetypes.sort_values(by='Frequency', ascending=False)
            clonetypes = clonetypes.reset_index()

            clonetypes['CloneId'] = [i for i in range(1, (clonetypes.shape[0]+1))]
            clonetypes['TRA_CDR3aa'] = clonetypes.TRA_CDR3aa.apply(lambda x: 'C'+str(x)+'F' if x != 'None' else 'None')
            clonetypes['TRB_CDR3aa'] = clonetypes.TRB_CDR3aa.apply(lambda x: 'C'+str(x)+'F' if x != 'None' else 'None')

            clonetypes = clonetypes.reindex(columns=list(['CloneId', 'TRA_CDR3aa', 'TRB_CDR3aa', 'Frequency', 'Proportion']))

            clonetypes.to_csv(f'{outdir}/clonetypes.tsv', sep='\t', index=None)

            vdj_sum_summary.append({
                'item': 'Estimated Number of Cells',
                'count': productive_cells_num,
                'total_count': np.nan,
            })

            loci = ['A', 'B']

            for locus in loci:
                tmp = int(clonetypes[clonetypes[f'TR{locus}_CDR3aa'] != 'None']['Frequency'].sum())

                vdj_sum_summary.append({
                    'item': f'Cells with TR{locus}',
                    'count': tmp,
                    'total_count': productive_cells_num,
                })

            paired_cell = int(clonetypes[(clonetypes['TRA_CDR3aa'] != 'None') & (clonetypes['TRB_CDR3aa'] != 'None')]['Frequency'].sum())

            vdj_sum_summary.append({
                'item': 'Cells with paired TRA and TRB',
                'count': paired_cell,
                'total_count': productive_cells_num,
            })
          

            for locus in loci:
                tmp = glob.glob(f'{ass_dir}/tracer/*/aligned_reads/*_TCR_{locus}.fastq')
                if len(tmp) != 0:
                    read_count = [get_umi_count(fq) for fq in tmp]
                    read_count.sort()
                    for i in range(len(read_count)):
                        if read_count[i] != 0:
                            idx = i
                            break
                    read_count = read_count[idx:]
                    median_tmp = int(np.median(read_count))
                    vdj_sum_summary.append({
                        'item': f'Median TR{locus} UMIs per cell',
                        'count': median_tmp,
                        'total_count': np.nan
                    })
                else:
                    vdj_sum_summary.append({
                    'item': f'Median TR{locus} UMIs per cell',
                    'count': 0,
                    'total_count': np.nan
                    })


        elif Seqtype == 'BCR':

            # barcode umi plot
            count_umi['mark'] = count_umi['cell_name'].apply(lambda x: "CB" if (x in productive_cells) else "UB")

            count_umi.to_csv(count_umi_file, sep='\t')        

            self.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_umi_file))


            # clone type table

            Proportion = []
            sum_c = clonetypes['Frequency'].sum()
            for f in list(clonetypes['Frequency']):
                p = f/sum_c
                p = p * 100
                p = round(p, 2)
                p = str(p) + '%'
                Proportion.append(p)
            clonetypes['Proportion'] = Proportion
            clonetypes = clonetypes.sort_values(by='Frequency', ascending=False)
            clonetypes = clonetypes.reset_index()

            clonetypes['CloneId'] = [i for i in range(1, (clonetypes.shape[0]+1))]
            clonetypes = clonetypes.reindex(columns=list(['CloneId', 'IGH_CDR3aa', 'IGL_CDR3aa', 'IGK_CDR3aa', 'Frequency', 'Proportion']))

            clonetypes.to_csv(f'{outdir}/clonetypes.tsv', sep='\t', index=None)


            vdj_sum_summary.append({
                    'item': 'Estimated Number of Cells',
                    'count': productive_cells_num,
                    'total_count': np.nan
            })

            loci = ['H', 'L', 'K']

            for locus in loci:
                tmp = int(clonetypes[clonetypes[f'IG{locus}_CDR3aa']!='None']['Frequency'].sum())

                vdj_sum_summary.append({
                        'item': f'Cells with IG{locus}',
                        'count': tmp,
                        'total_count': productive_cells_num
                })

            paired_H_L = int(clonetypes[(clonetypes['IGH_CDR3aa']!='None') & (clonetypes['IGL_CDR3aa']!='None')]['Frequency'].sum())

            vdj_sum_summary.append({
                    'item': 'Cells with paired IGH and IGL',
                    'count': paired_H_L,
                    'total_count': productive_cells_num
            })

            paired_H_K = int(clonetypes[(clonetypes['IGH_CDR3aa']!='None') & (clonetypes['IGK_CDR3aa']!='None')]['Frequency'].sum())

            vdj_sum_summary.append({
                    'item': 'Cells with paired IGH and IGK',
                    'count': paired_H_K,
                    'total_count': productive_cells_num
            })


            for locus in loci:
                tmp = glob.glob(f'{ass_dir}/bracer/*/aligned_reads/*_BCR_{locus}.fastq')
                if len(tmp) != 0:
                    read_count = [get_umi_count(fq) for fq in tmp]    
                    read_count.sort()
                    for i in range(len(read_count)):
                        if read_count[i] != 0:
                            idx = i
                            break
                    read_count = read_count[idx:]
                    median_tmp = int(np.median(read_count))
                    vdj_sum_summary.append({
                    'item': f'Median IG{locus} UMIs per cell',
                    'count': median_tmp,
                    'total_count': np.nan
                    })
                else:
                    vdj_sum_summary.append({
                    'item': f'Median IG{locus} UMIs per cell',
                    'count': 0,
                    'total_count': np.nan
                    })

        df = pd.DataFrame(vdj_sum_summary, 
            columns=['item', 'count', 'total_count'])

        utils.gen_stat(df, stat_file)

        # clonetype table

        title = 'Clonetypes'
        table_dict = self.get_table(title, 'clonetypes_table', clonetypes)

        self.add_data_item(table_dict=table_dict)

        self.clean_up()

@utils.add_log
def vdj_sum(args):
    step_name = 'vdj_sum'
    vdj_sum_obj = Vdj_sum(args, step_name)
    vdj_sum_obj.run()


def get_opts_vdj_sum(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--ass_dir', help='assemble dir', required=True)
        parser.add_argument('--fastq_dir', help='dir contains fastq', required=True)
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)



