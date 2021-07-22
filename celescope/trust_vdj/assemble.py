import os
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from celescope.tools import utils
from celescope.tools.cellranger3 import get_plot_elements
from celescope.tools.step import Step, s_common

SOFTWARE = '/SGRNJ03/randd/zhouxin/software/TRUST4/'

class Assemble(Step):
    """
    Features

    - Assemble TCR/BCR seq data.

    Output

    - `05.assemble/{sample}_toassemble.fq` Reads to assemble.
    - `05.assemble/{sample}_toassemble_bc.fa` Barcodes to assemble.
    - `05.assemble/{sample}_cdr3.out` All assembled CDR3 output.
    - `05.assemble/{sample}_barcode_report.tsv` Record chain information in each barcode.
    - `05.assemble/{sample}_annot.fa` Assembled annotated contig sequences.
    - `05.assemble/{sample}_assembled_reads.fa` Assembled raw reads.
    - `05.assemble/{sample}_report.tsv` Record assembled CDR3 types and count.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.sample = args.sample
        self.species = args.species
        self.speed_up = args.speed_up
        self.Seqtype = args.Seqtype
        self.cells = args.cells
        self.cb_stat = args.cb_stat
        
        if self.Seqtype == 'TCR':
            self.string = 't'
            self.chain = ['TRA', 'TRB']
            self.paired_groups = ['TRA_TRB']
        elif self.Seqtype == 'BCR':
            self.string = 'b'
            self.chain = ['IGH', 'IGL', 'IGK']
            self.paired_groups = ['IGH_IGL', 'IGH_IGK']
            
        # res out prefix 
        self.final_out = f'{self.outdir}/outs'
        
        # trust out file
        self.report = f'{self.outdir}/{self.sample}_barcode_report.tsv'
        self.cdr3out = f'{self.outdir}/{self.sample}_cdr3.out'
        self.fa = f'{self.outdir}/{self.sample}_annot.fa'
        self.assembled_reads = f'{self.outdir}/{self.sample}_assembled_reads.fa'
        self.assign_file = f'{self.outdir}/{self.sample}_assign.out'
        self.count_file = f'{self.outdir}/count.txt'

        # format out
        if not os.path.exists(self.final_out):
            os.makedirs(self.final_out)
        
        self.all_out_rep_prefix = f'{self.final_out}/all'
        self.filter_out_rep_prefix  = f'{self.final_out}/filter'
        self.all_tmp_rep = f'{self.all_out_rep_prefix}_{self.string}.csv'
        self.filter_tmp_rep = f'{self.filter_out_rep_prefix}_{self.string}.csv'
        self.all_rep = f'{self.all_out_rep_prefix}_contig_annotations.csv'
        self.filter_rep = f'{self.filter_out_rep_prefix}ed_contig_annotations.csv'
        
        # inde
        species = self.species

        self.index_file = f'{SOFTWARE}/index/{species}/{species}_ref.fa'
        self.ref = f'{SOFTWARE}/index/{species}/{species}_IMGT+C.fa'
        
        # summary
        self.assemble_summary = []
        self.cells_num = 0
        self.mapping_reads = 0
        
        
        # match_dir
        self.match_bool = True
        if (not args.match_dir) or (args.match_dir == "None"):
            self.match_bool = False
        if self.match_bool:
            self.match_cell_barcodes, _match_cell_number = utils.read_barcode_file(
                args.match_dir)
        
        
    @utils.add_log
    def mapping(self):
        cmd = (
            f'{SOFTWARE}/fastq-extractor -t {self.thread} -f {SOFTWARE}/index/{self.species}/{self.species}_ref.fa '
            f'-o {self.outdir}/{self.sample}_raw '  
            f'--barcodeStart 0 --barcodeEnd 23 --umiStart 24 --umiEnd 31 '
            f'-u {self.fq2} --barcode {self.fq1} --UMI {self.fq1}'
        )
        Assemble.mapping.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

        
    @utils.add_log
    def cutoff(self):
        dic = defaultdict(set)
        read_dic = defaultdict(list)
        with pysam.FastxFile(f'{self.outdir}/{self.sample}_raw.fq', 'r') as fa:
            for entry in fa:
                name = entry.name
                attrs = name.split('_')
                cb = attrs[0]
                umi = attrs[1]
                dic[cb].add(umi)
                read_dic[cb].append(entry)
                self.mapping_reads+=1
                
        df = pd.DataFrame()
        df['barcode'] = list(dic.keys())
        df['UMI'] = [len(dic[i]) for i in list(dic.keys())]
        df_sort = df.sort_values(by='UMI', ascending=False)
        #df_sort = df_sort.reset_index(drop=True)
        df_sort.to_csv(self.count_file, sep='\t', index=False)
        
        UMI_num = int(self.cells)
        rank = int(UMI_num / 100)
        rank_UMI = df_sort.iloc[rank, :]['UMI']
        UMI_min = int(rank_UMI / 10)

        df_umi_filtered = df_sort[df_sort['UMI'] >= UMI_min]
        filtered_barcodes = df_umi_filtered['barcode'].tolist()
        
        if self.match_bool:
            barcodes = set(filtered_barcodes).intersection(set(self.match_cell_barcodes))
        else:
            barcodes = filtered_barcodes
        Assemble.cutoff.logger.info(f'get {len(barcodes)} cells to assemble')
        
        new_fq = open(f'{self.outdir}/{self.sample}_toassemble.fq', 'w')
        new_cb = open(f'{self.outdir}/{self.sample}_toassemble_bc.fa', 'w')
        new_umi = open(f'{self.outdir}/{self.sample}_toassemble_umi.fa', 'w')
        read_count = 0
        barcode_count = set()
        for barcode in barcodes:
            read_list = read_dic[barcode]
            for entry in read_list:
                read_id = entry.name
                attrs = read_id.split('_')
                cb = attrs[0]
                umi = attrs[1]
                new_fq.write(f'{str(entry)}\n')
                new_cb.write(f'>{read_id}\n{cb}\n')
                new_umi.write(f'>{read_id}\n{umi}\n')
                read_count += 1
                barcode_count.add(cb)
                if read_count % 1000000 == 0:
                    Assemble.cutoff.logger.info(f'processed {read_count} reads')
        Assemble.cutoff.logger.info(f'confirmed {len(barcode_count)} cells to assemble')
            
        new_fq.close()
        new_cb.close()
        new_umi.close()
            
              
    @utils.add_log
    def run_assemble(self):

        string1 = ''
        if self.speed_up:
            string1 = '--repseq '
        cmd = (
            f'source activate zhouxinT; '
            f'{SOFTWARE}/run-trust4 -t {self.thread} '
            f'-u {self.fq2} '
            f'--barcode {self.fq1} '
            f'--barcodeRange 0 23 + '
            f'--UMI {self.fq1} '
            f'--umiRange 24 31 + '
            f'-f {self.index_file} '
            f'--ref {self.ref} '
            f'{string1}'
            f'-o {self.sample} --od {self.outdir} --outputReadAssignment '
            f'--stage 1'
        )
        Assemble.run_assemble.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
        
    @utils.add_log
    def get_all_rep(self):
        cmd = (
           f'perl {SOFTWARE}/trust-barcoderep-to-10X.pl '
           f'{self.report} '
           f'{self.all_out_rep_prefix} '
        )
        subprocess.check_call(cmd, shell=True)
        Assemble.get_all_rep.logger.info(cmd)
        os.system(f'mv {self.all_tmp_rep} {self.all_rep}')


    @utils.add_log
    def get_len(self):
        all_rep = pd.read_csv(self.all_rep, sep=',')
        with pysam.FastaFile(self.fa) as fh:
            all_rep['length'] = all_rep['contig_id'].apply(lambda x: len(fh.fetch(x)))
            
            all_rep[['reads', 'umis']] = all_rep[['reads', 'umis']].astype(int)

            all_rep.to_csv(self.all_rep, sep=',', index=False)
            
                 
    @utils.add_log
    def get_filter_rep(self):
        all_rep = pd.read_csv(self.all_rep, sep=',')
        filtered_rep = all_rep[(all_rep['productive']==True) & (all_rep['full_length']==True)]
        filtered_rep.to_csv(self.filter_rep, sep=',', index=False)
        
        # get stat
        cell_barcodes = filtered_rep['barcode'].tolist()
        cells_num = len(set(cell_barcodes))
        self.cells_num = cells_num
        self.assemble_summary.append({
            'item': 'Estimated Number of Cells',
            'count': cells_num,
            'total_count': np.nan
        })
        
        df_umi = pd.read_csv(self.count_file, sep='\t')
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in cell_barcodes else 'UB')
        df_umi.to_csv(self.count_file, sep='\t', index=False)
        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(self.count_file))
           
        
    @utils.add_log
    def get_clonetypes(self):
        data = pd.read_csv(self.filter_rep, sep=',')
        with open(f'{self.final_out}/clonetypes.csv', 'w') as fh:
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

        df = pd.read_csv(f'{self.final_out}/clonetypes.csv', sep=',')
        df = df.groupby(['cdr3s_aa', 'cdr3s_nt'], as_index=False).agg({'barcode': 'count'})
        df = df.rename(columns={'barcode': 'frequency'})
        df = df.sort_values(by='frequency', ascending=False)
        sum_c = df['frequency'].sum()
        df['proportion'] = df['frequency'].apply(lambda x: x/sum_c)
        df['clonotype_id'] = [f'clonetype{i}' for i in range(1, df.shape[0]+1)]

        df = df.reindex(columns=['clonotype_id', 'frequency', 'proportion', 'cdr3s_aa', 'cdr3s_nt'])
        df.to_csv(f'{self.final_out}/clonetypes.csv', sep=',', index=False)
        
    @utils.add_log
    def gen_table(self):
        df = pd.read_csv(f'{self.final_out}/clonetypes.csv', sep=',')
        df_table = pd.DataFrame()
        def function_x(x):
            l = x.split(';')
            if len(l)==1:
                if l[0].startswith('TRA'):
                    return l[0], 'None'
                elif l[0].startswith('TRB'):
                    return 'None', l[0]
            elif len(l)==2:
                return l[0], l[1]
        df_table['Clonotype ID'] = df['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        df_table['Chain1'] = df['cdr3s_aa'].apply(lambda x: f'{function_x(x)[0]}')
        df_table['Chain2'] = df['cdr3s_aa'].apply(lambda x: f'{function_x(x)[1]}')
        df_table['Frequency'] = df['frequency']
        df_table['Proportion'] = df['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        title = 'Clonetypes'

        table_dict = self.get_table(title, 'clonetypes_table', df_table)

        self.add_data_item(table_dict=table_dict)
        
    @utils.add_log
    def get_full_len_assembly(self):
        cmd = (
            f'perl {SOFTWARE}/GetFullLengthAssembly.pl '
            f'{self.fa} '
            f'> {self.final_out}/full_len_contig.fa '
        )
        subprocess.check_call(cmd, shell=True)
            
    @utils.add_log
    def get_stat_file(self):
        with open(self.cb_stat, 'r') as fh:
            for line in fh:
                if 'Raw Reads:' in line:
                    count = int(line.split(':')[1].replace(',',''))
                    break

        self.assemble_summary.append({
            'item': 'Mean Read Pairs per Cell',
            'count': int(count/self.cells_num),
            'total_count': np.nan
        })
        
        fa = pysam.FastaFile(self.assembled_reads)
        assembled_reads = fa.nreferences
        self.assemble_summary.append({
            'item': 'Mean Used Read Pairs per Cell',
            'count': int(assembled_reads/self.cells_num),
            'total_count': np.nan
        })
        self.assemble_summary.append({
            'item': 'Fraction Reads in Cells',
            'count': assembled_reads,
            'total_count': count
        })
        
        data = pd.read_csv(self.filter_rep, sep=',')
        data = data[['barcode', 'chain', 'contig_id']]
        assign_reads = pd.read_csv(self.assign_file, sep='\t', header=None)
        assign_reads = assign_reads.rename(columns={0: 'read_id', 1: 'contig_id'})
        assign_df = pd.merge(data, assign_reads, on='contig_id', how='left').fillna('None')
        self.assemble_summary.append({
            'item': 'Reads Mapped to Any V(D)J Gene',
            'count': self.mapping_reads,
            'total_count': count
        })
        for c in self.chain:
            # self.assemble_summary.append({
            #     'item': f'Reads Mapped to {c}', 
            #     'count': df_c.shape[0], 
            #     'total_count': count
            # })
            dic = defaultdict(set)
            tmp = assign_df[assign_df['chain']==c]
            contig_ids = data[data['chain']==c]['contig_id'].tolist()
            for i in contig_ids:
                cell_bc = i.split("_")[0]
                tmp_ = tmp[tmp['contig_id']==i]
                for read_id in tmp_['read_id'].tolist():
                    attrs = read_id.split('_')
                    if len(attrs)>1:
                        umi = attrs[1]
                        dic[cell_bc].add(umi)
            umi_count = [len(dic[i]) for i in list(dic.keys())]
            mid = int(np.median(umi_count))
            item = f'Median {c} UMIs per cell'
            self.assemble_summary.append({
                'item': item,
                'count': mid,
                'total_count': np.nan
            })
        
        assemble_summary = self.assemble_summary
        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(assemble_summary, columns=['item', 'count', 'total_count'])

        utils.gen_stat(sum_df, stat_file)
        
    def run(self):
        #if not os.path.exists(f'{self.outdir}/{self.sample}_raw.fq'):
        self.mapping()
            
        self.cutoff()

        if not os.path.exists(f'{self.outdir}/{self.sample}_annot.fa'):
            self.run_assemble()

        self.get_all_rep()
        self.get_len()
        self.get_filter_rep()
        self.get_clonetypes()
        self.gen_table()
        self.get_full_len_assembly()
        self.get_stat_file()
        self.clean_up()

def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='R1 reads from match step', required=True)
        parser.add_argument('--fq2', help='R2 reads from match step', required=True)
        parser.add_argument('--cb_stat', help='barcode stat file', required=True)
        parser.add_argument('--match_dir', help='scRNA-seq results directory', default=None)
    parser.add_argument('--cells', help='expected cells num', default=3000)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)
    parser.add_argument('--speed_up', help='speed assemble for TCR/BCR seq data', action='store_true')
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)      








