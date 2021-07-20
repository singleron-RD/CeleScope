import os
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.cellranger3 import get_plot_elements
from celescope.tools.step import Step, s_common
from celescope.vdj10X.__init__ import ref_dict, soft_dict

pd.set_option('max_colwidth',200)

def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())

def get_vj_annot(df, chains, pairs):
    fl_pro_pair_df = pd.DataFrame(df[(df['full_length']==True)&(df['productive']==True)].drop_duplicates(['barcode', 'chain']).barcode.value_counts())
    fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']==2]
    l = []
    cell_nums = len(set(df['barcode'].tolist()))
    l.append({
        'item': 'Cells With Productive V-J Spanning Pair',
        'count': fl_pro_pair_df.shape[0],
        'total_count': cell_nums
    })
    for p in pairs:
        chain1 = p.split('_')[0]
        chain2 = p.split('_')[1]
        cbs1 = set(df[(df['productive']==True)&(df['chain']==chain1)].barcode.tolist())
        cbs2 = set(df[(df['productive']==True)&(df['chain']==chain2)].barcode.tolist())
        paired_cbs = len(cbs1.intersection(cbs2))
        l.append({
            'item': f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
            'count': paired_cbs,
            'total_count': cell_nums
        })
    for c in chains:
        l.append({
            'item': f'Cells With {c} Contig',
            'count': len(set(df[df['chain']==c].barcode.tolist())),
            'total_count': cell_nums
        })
        l.append({
            'item': f'Cells With CDR3-annotated {c} Contig',
            'count': len(set(df[(df['chain']==c)&(df['cdr3']!='None')].barcode.tolist())),
            'total_count': cell_nums
        })
        l.append({
            'item': f'Cells With V-J Spanning {c} Contig',
            'count': len(set(df[(df['full_length']==True)&(df['chain']==c)].barcode.tolist())),
            'total_count': cell_nums
        })
        l.append({
            'item': f'Cells With Productive {c} Contig',
            'count': len(set(df[(df['productive']==True)&(df['chain']==c)].barcode.tolist())),
            'total_count': cell_nums
        })
    # l.append({
    #     'item': 'Cells With TRA Contig',
    #     'count': len(set(df[df['chain']=='TRA'].barcode.tolist())),
    #     'total_count': cell_nums
    # })
    # l.append({
    #     'item': 'Cells With TRB Contig',
    #     'count': len(set(df[df['chain']=='TRB'].barcode.tolist())),
    #     'total_count': cell_nums
    # })
    # l.append({
    #     'item': 'Cells With CDR3-annotated TRA Contig',
    #     'count': len(set(df[(df['cdr3']!='None')&(df['chain']=='TRA')].barcode.tolist())),
    #     'total_count': cell_nums
    # })
    # l.append({
    #     'item': 'Cells With CDR3-annotated TRB Contig',
    #     'count': len(set(df[(df['cdr3']!='None')&(df['chain']=='TRB')].barcode.tolist())),
    #     'total_count': cell_nums
    # })
    # l.append({
    #     'item': 'Cells With V-J Spanning TRA Contig',
    #     'count': len(set(df[(df['full_length']==True)&(df['chain']=='TRA')].barcode.tolist())),
    #     'total_count': cell_nums
    # })
    # l.append({
    #     'item': 'Cells With V-J Spanning TRB Contig',
    #     'count': len(set(df[(df['full_length']==True)&(df['chain']=='TRB')].barcode.tolist())),
    #     'total_count': cell_nums
    # })
    # l.append({
    #     'item': 'Cells With Productive TRA Contig',
    #     'count': len(set(df[(df['productive']==True)&(df['chain']=='TRA')].barcode.tolist())),
    #     'total_count': cell_nums
    # })
    # l.append({
    #     'item': 'Cells With Productive TRB Contig',
    #     'count': len(set(df[(df['productive']==True)&(df['chain']=='TRB')].barcode.tolist())),
    #     'total_count': cell_nums
    # })
    return l
     
class Assemble(Step):
    
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        
        # common parameters
        self.mem = args.mem
        self.species = args.species
        self.soft = args.soft
        self.thread = args.thread
        self.fqs_dir = args.fqs_dir
        self.match_dir = args.match_dir
        self.seqtype = args.seqtype
        
        # chain parameters
        if self.seqtype == 'TCR':
            self.chains = ['TRA', 'TRB']
            self.pair = ['TRA_TRB']
        elif self.seqtype == 'BCR':
            self.chains = ['IGH', 'IGL', 'IGK']
            self.pair = ['IGH_IGL', 'IGH_IGK']
        
        # input
        self.barcode_dic = args.barcode_dic
        self.all_bam = f'{self.outdir}/{self.sample}/outs/all_contig.bam'
        self.filter_contig = f'{self.outdir}/{self.sample}/outs/filtered_contig_annotations.csv'
        self.filter_fa = f'{self.outdir}/{self.sample}/outs/filtered_contig.fasta'
        
        # output
        self.all_outs = f'{self.outdir}/all'
        self.out_fa = f'{self.all_outs}/filtered_contig.fasta'
        self.match_outs = f'{self.outdir}/match'
        self.match_fa = f'{self.match_outs}/match_contig.fasta'
        
        if not os.path.exists(self.all_outs):
            os.mkdir(self.all_outs)
        if not os.path.exists(self.match_outs):
            os.mkdir(self.match_outs)
        
        #other parameters
        self.cwd_path = os.getcwd()
        self.outs = f'{self.outdir}/{self.sample}/outs'
        
        # get match barcodes from match_dir
        self.match_bool = True
        if (not args.match_dir) or (args.match_dir == "None"):
            self.match_bool = False
        if self.match_bool:
            self.match_cell_barcodes, _match_cell_number = utils.read_barcode_file(
                args.match_dir)
             
    @utils.add_log
    def run_assemble(self):
        ref_path = ref_dict[self.soft][self.species]
        soft_path = soft_dict[self.soft]
        cmd = (
            f'{soft_path} vdj '
            f'--id={self.sample} '
            f'--reference={ref_path} '
            f'--fastqs={self.cwd_path}/{self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )
        Assemble.run_assemble.logger.info(cmd)
        with open(f'{self.outdir}/{self.sample}_vdj_10X.sh', 'w') as f:
            f.write(cmd)
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
    
    @utils.add_log   
    def gen_result(self):
        # set work directory
        os.chdir(self.cwd_path)
        
        # get barcode correspond dictionary
        barcode_df = pd.read_csv(self.barcode_dic, sep='\t', index_col=1)
        barcode_dict = barcode_df.to_dict()['sgr']
        
        # get sum_metrics
        sum_dict = pd.read_csv(f'{self.outs}/metrics_summary.csv', sep=',', index_col=None)
        sum_dict = sum_dict.T.to_dict()
        read_count = int(sum_dict[0]["Number of Read Pairs"].replace(',', ''))
        
        # gen clone table
        raw_clonotypes = pd.read_csv(f'{self.outs}/clonotypes.csv', sep=',', index_col=None)
        raw_clonotypes['ClonotypeID'] = raw_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        raw_clonotypes['Frequency'] = raw_clonotypes['frequency']
        raw_clonotypes['Proportion'] = raw_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        raw_clonotypes['CDR3_aa'] = raw_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))
        
        # write table
        title = 'Clonetypes'
        table_dict = self.get_table(title, 'clonetypes_table', raw_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']])
        self.add_data_item(table_dict=table_dict)
        all_clonotypes = raw_clonotypes[['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion']]
        all_clonotypes.to_csv(f'{self.all_outs}/clonotypes.csv', sep=',', index=False)

        # convert barcode to sgr_cbs
        ## filtered_contig_annotations.csv
        filter_contig = pd.read_csv(f'{self.filter_contig}', sep=',', index_col=None)
        filter_contig['barcode'] = filter_contig['barcode'].apply(lambda x: reversed_compl(barcode_dict[x.split('-')[0]]))
        filter_contig['contig_id'] = filter_contig['contig_id'].apply(lambda x: reversed_compl(barcode_dict[x.split('-')[0]])+'_'+x.split('_')[1]+'_'+x.split('_')[2])
        filter_contig.to_csv(f'{self.all_outs}/filtered_contig_annotations.csv', sep=',', index=False)
        ## filtered_contig.fasta
        in_fa = pysam.FastxFile(self.filter_fa)
        out_fa = open(self.out_fa, 'w')
        for entry in in_fa:
            name = entry.name
            seq = entry.sequence
            attrs = name.split('_')
            new_name = reversed_compl(barcode_dict[attrs[0].split('-')[0]]) + '_' + attrs[1] + '_' + attrs[2]
            out_fa.write(f'>{new_name}\n{seq}\n')
        out_fa.close()
             
        # get common stat file
        common_summary = []
        cell_nums = len(set(filter_contig['barcode'].tolist()))
        common_summary.append({
            'item': 'Estimated Number of Cells',
            'count': cell_nums,
            'total_count': np.nan
        })
        common_summary.append({
            'item': 'Reads Mapped to Any V(D)J Gene',
            'count': int(read_count * (float(sum_dict[0]['Reads Mapped to Any V(D)J Gene'].strip('%'))/100)), 
            'total_count': read_count
        })
        for c in self.chains:
            common_summary.append({
                'item': f'Reads Mapped to {c}', 
                'count': int(read_count * (float(sum_dict[0][f'Reads Mapped to {c}'].strip('%'))/100)), 
                'total_count': read_count
            })
        # common_summary.append({
        #     'item': 'Reads Mapped to TRA',
        #     'count': int(read_count * (float(sum_dict[0]['Reads Mapped to TRA'].strip('%'))/100)),
        #     'total_count': read_count,
        # })
        # common_summary.append({
        #     'item': 'Reads Mapped to TRB',
        #     'count': int(read_count * (float(sum_dict[0]['Reads Mapped to TRB'].strip('%'))/100)),
        #     'total_count': read_count
        # })
        common_summary.append({
            'item': 'Fraction Reads in Cells',
            'count': int(read_count * (float(sum_dict[0]['Fraction Reads in Cells'].strip('%'))/100)),
            'total_count': read_count
        })
        for c in self.chains:
            common_summary.append({
                'item': f'Median used {c} UMIs per Cell',
                'count': int(filter_contig[filter_contig['chain']==c]['umis'].median()),
                'total_count': np.nan
            })
        # common_summary.append({
        #     'item': 'Median used TRA UMIs per Cell',
        #     'count': int(filter_contig[filter_contig['chain']=='TRA']['umis'].median()), 
        #     'total_count': np.nan
        # })
        # common_summary.append({
        #     'item': 'Median used TRB UMIs per Cell',
        #     'count': int(filter_contig[filter_contig['chain']=='TRB']['umis'].median()), 
        #     'total_count': np.nan
        # })
        # get all results stat_file
        all_summary = get_vj_annot(filter_contig, self.chains, self.pair)   

        # count umi and plot
        all_bam = pysam.AlignmentFile(self.all_bam)
        dic_umi = defaultdict(set)
        
        for read in all_bam:
            cb = read.get_tag('CB')
            umi = read.get_tag('UB')
            new_cb = reversed_compl(barcode_dict[cb.split('-')[0]])
            dic_umi[new_cb].add(umi)
            
        df_umi = pd.DataFrame()
        df_umi['barcode'] = list(dic_umi.keys())
        df_umi['UMI'] = [len(dic_umi[i]) for i in dic_umi]
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        sgr_cbs = set(filter_contig['barcode'].tolist())
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in sgr_cbs else 'UB')
        df_umi.to_csv(f'{self.all_outs}/count.txt', sep='\t', index=False)
        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(f'{self.all_outs}/count.txt'))
        
        # match
        df_sgr = pd.DataFrame(self.match_cell_barcodes, columns=['barcode'])
        df_match = pd.merge(df_sgr, filter_contig, on='barcode', how='inner')
        
        # get match summary
        match_summary = get_vj_annot(df_match, self.chains, self.pair)
        match_summary.insert(0, {
            'item': 'Cells match with scRNA-seq analysis',
            'count': len(set(df_match['barcode'].tolist())),
            'total_count': np.nan
        })
        
        # gen match results
        # df_match['contig_id'] = df_match['contig_id'].apply(lambda x: reversed_compl(barcode_dict[x.split('-')[0]]))
        df_match.to_csv(f'{self.match_outs}/match_contigs.csv', sep=',', index=False)
        
        # get match clonotypes
        df_match = df_match[df_match['productive']==True]
        df_match['chain_cdr3aa'] = df_match[['chain', 'cdr3']].apply(':'.join, axis=1)
        match_cbs = set(df_match['barcode'].tolist())
        match_clonotypes = open(f'{self.match_outs}/match_clonotypes.csv', 'w')
        match_clonotypes.write('barcode\tcdr3s_aa\n')
        for cb in match_cbs:
            temp = df_match[df_match['barcode']==cb]
            temp = temp.sort_values(by='chain', ascending=True)
            chain_list = temp['chain_cdr3aa'].tolist()
            chain_str = ';'.join(chain_list)
            # print(chain_str)
            match_clonotypes.write(f'{cb}\t{chain_str}\n')
        match_clonotypes.close()
            
        df_match_clonetypes = pd.read_csv(f'{self.match_outs}/match_clonotypes.csv', sep='\t', index_col=None)
        df_match_clonetypes = df_match_clonetypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
        df_match_clonetypes = df_match_clonetypes.rename(columns={'barcode': 'frequency'})
        sum_f = df_match_clonetypes['frequency'].sum()
        df_match_clonetypes['proportion'] = df_match_clonetypes['frequency'].apply(lambda x: x/sum_f)
        df_match_clonetypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_match_clonetypes.shape[0]+1)]
        df_match_clonetypes = df_match_clonetypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
        df_match_clonetypes = df_match_clonetypes.sort_values(by='frequency', ascending=False)
        df_match_clonetypes.to_csv(f'{self.match_outs}/match_clonotypes.csv', sep=',', index=False)
        
        fa = pysam.FastxFile(self.filter_fa)
        match_fa = open(self.match_fa, 'w')
        for entry in fa:
            name = entry.name
            attrs = name.split('_')
            cb = attrs[0].split('-')[0]
            new_cb = reversed_compl(barcode_dict[cb])
            if new_cb in self.match_cell_barcodes:
                new_name = new_cb + '_' + attrs[1] + '_' + attrs[2]
                seq = entry.sequence
                match_fa.write(f'>{new_name}\n{seq}\n')
        
        #get assemble_summary
        assemble_summary = common_summary + all_summary + match_summary
        
        # get stat file
        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(assemble_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)
        
    def run(self):
        
        if not os.path.exists(self.outs):
            self.run_assemble()
        self.gen_result()
        self.clean_up()
            
def assemble(args):
    step_name = f'{args.seqtype}_assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    parser.add_argument('--species', help='species', choices=['hs','mmu'], required=True)
    parser.add_argument('--soft', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'], 
        default='4.0.0')
    parser.add_argument('--mem', help='memory (G)', default=10)
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir', required=True)
        parser.add_argument('--barcode_dic', help='10X barcode correspond sgr barcode', required=True)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
    return parser
