from celescope.tools import utils
from celescope.tools.step import Step, s_common
import subprocess
import os
import pysam
from collections import defaultdict
import pandas as pd
from celescope.tools.cellranger3 import get_plot_elements
import numpy as np

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
        self.all_rep = f'{self.all_out_rep_prefix}ed_contig_annotations.csv'
        self.filter_rep = f'{self.filter_out_rep_prefix}ed_contig_annotations.csv'
        
        # inde
        species = self.species

        self.index_file = f'{SOFTWARE}/index/{species}/{species}_ref.fa'
        self.ref = f'{SOFTWARE}/index/{species}/{species}_IMGT+C.fa'
        
        # summary
        self.assemble_summary = []
        self.cells_num = 0
        self.contig_df = pd.DataFrame()
        
        
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
        barcodes = df_umi_filtered['barcode'].tolist()
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
            
        subprocess.check_call(f'rm {self.outdir}/{self.sample}_raw.fq', shell=True)
        subprocess.check_call(f'rm {self.outdir}/{self.sample}_raw_bc.fa', shell=True)
        subprocess.check_call(f'rm {self.outdir}/{self.sample}_raw_umi.fa', shell=True)
            
              
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
        os.system(f'rm {self.final_out}/all_*.csv')
        
    @utils.add_log
    def get_filter_rep(self):
        out_rep = f'{self.outdir}/{self.sample}_filter_report.tsv'
        cmd = (
            f'perl {SOFTWARE}/trust-barcoderep.pl '
            f'{self.cdr3out} '
            f'-a {self.fa} '
            f'--noImputation > '
            f'{out_rep}; '
            f'perl {SOFTWARE}/trust-barcoderep-to-10X.pl '
            f'{out_rep} '
            f'{self.filter_out_rep_prefix} '
        )
        subprocess.check_call(cmd, shell=True)
        Assemble.get_filter_rep.logger.info(cmd)
        os.system(f'mv {self.filter_tmp_rep} {self.filter_rep}')
        os.system(f'rm {self.final_out}/filter_*.csv')

        os.remove(out_rep)
        
    @utils.add_log
    def get_len(self):
        with pysam.FastxFile(self.fa, 'r') as fh:
            res = defaultdict(list)
            for entry in fh:
                name = entry.name
                comment = entry.comment
                attrs = comment.split(' ')
                V_gene = attrs[2]
                D_gene = attrs[3]
                J_gene = attrs[4]
                if V_gene != '*':
                    chain = V_gene[:3]
                elif J_gene != '*':
                    chain = J_gene[:3]
                elif D_gene != '*':
                    chain = D_gene[:3]
                else:
                    chain = 'None'
                length = len(entry.sequence)
                res['contig_id'].append(name)
                res['length'].append(length)
                res['chain'].append(chain)
            
            df = pd.DataFrame(res, columns=list(res.keys()))
            self.contig_df = df
            df = df.set_index(['contig_id'])
            
            filter_res = pd.read_csv(self.filter_rep, sep=',')
            filter_res['length'] = filter_res['contig_id'].apply(lambda x: df.loc[x, 'length'])
            filter_res[['reads', 'umis']] = filter_res[['reads', 'umis']].astype(int)

            filter_res.to_csv(self.filter_rep, sep=',', index=False)
            
            all_res = pd.read_csv(self.all_rep, sep=',')
            all_res['length'] = all_res['contig_id'].apply(lambda x: df.loc[x, 'length'])
            all_res[['reads', 'umis']] = all_res[['reads', 'umis']].astype(int)

            all_res.to_csv(self.all_rep, sep=',', index=False)
            
            cell_barcodes = filter_res['barcode'].tolist()
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
        data = pd.read_csv(f'{self.filter_rep}', sep=',')
        data = data[data['productive']==True]
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
                
            fh.close()

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
    def get_full_len_assembly(self):
        cmd = (
            f'perl {SOFTWARE}/GetFullLengthAssembly.pl '
            f'{self.fa} '
            f'> {self.final_out}/full_len_contig.fa '
        )
        subprocess.check_call(cmd, shell=True)
        
        df = pd.read_csv(f'{self.filter_rep}', sep=',')
        df_pro = df[df['productive']==True]
        contig_ids = df_pro['contig_id'].tolist()
        
        with pysam.FastxFile(f'{self.final_out}/full_len_contig.fa', 'r') as full_len_fa:
            pro_full_len_fa = open(f'{self.final_out}/productive_full_len_contig.fa', 'w')
            for entry in full_len_fa:
                name = entry.name
                if name in contig_ids:
                    pro_full_len_fa.write(str(entry) + '\n')
                    
            pro_full_len_fa.close()
            
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
        assign_df = pd.merge(self.contig_df, assign_reads, on='contig_id', how='outer')
        self.assemble_summary.append({
            'item': 'Reads Mapped to Any V(D)J Gene',
            'count': assign_df.shape[0],
            'total_count': count
        })
        for c in self.chain:
            df_c = assign_df[assign_df['chain']==c]
            self.assemble_summary.append({
                'item': f'Reads Mapped to {c}', 
                'count': df_c.shape[0], 
                'total_count': count
            })
            dic = defaultdict(set)
            tmp = data[data['chain']==c]
            contig_ids = tmp['contig_id'].tolist()
            for i in contig_ids:
                tm = assign_reads[assign_reads['contig_id']==i]
                barcode_ = i.split('_')[0]
                for read_id in tm['read_id'].tolist():
                    attrs = read_id.split('_')
                    cb = attrs[0]
                    umi = attrs[1]
                    if cb == barcode_:
                        dic[barcode_].add(umi)
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
        self.mapping()
        self.cutoff()
        self.run_assemble()
        self.get_all_rep()
        self.get_filter_rep()
        self.get_len()
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
    parser.add_argument('--cells', help='expected cells num', default=3000)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)
    parser.add_argument('--speed_up', help='speed assemble for TCR/BCR seq data', action='store_true')
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)      








