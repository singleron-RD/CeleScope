import pysam
from collections import defaultdict
import os
import argparse
import datetime
import pandas as pd
from Bio.Seq import Seq
import glob
import re
import numpy as np
from celescope.tools import utils
from celescope.tools.Step import Step, s_common
import glob
from celescope.tools.cellranger3 import get_plot_elements
import json


def tpm_count(ass_dir):
	rec = pd.read_csv(f'{ass_dir}/tracer/filtered_TCRAB_summary/recombinants.txt', sep='\t')  # ass_dir outdir/sample/04.go_assemble
	productive = rec[rec['productive'] == True]
	productive['TPM'] = ''
	indx = list(productive.index)
	for i in indx:
		cell_name = productive.at[i, 'cell_name']
		rec_id = productive.at[i, 'recombinant_id']
		with open(f'{ass_dir}/tracer/{cell_name}/expression_quantification/abundance.tsv') as tsvf:
			for line in tsvf:
				if rec_id in line:
					line = line.rstrip()
					line = line.split('\t')
					tpm = float(line[4])
					productive.loc[i, 'TPM'] = tpm
	
	return productive


def filtering(type, ass_dir, outdir):

	if not os.path.exists(outdir):
		os.makedirs(outdir)

	if type == 'TCR':
		data = tpm_count(ass_dir)
		cell_name = set(list(data['cell_name']))
		filtered = pd.DataFrame()
		for name in cell_name:
			count_data = data[data['cell_name'] == name]
			tra = count_data[count_data['locus'] == 'A']
			trb = count_data[count_data['locus'] == 'B']
			if tra.empty is not True:
				tra = tra.sort_values(by='TPM', ascending=False)
				tra = tra.head(1)
				filtered = filtered.append(tra, ignore_index=True)
			if trb.empty is not True:
				trb = trb.sort_values(by='TPM', ascending=False)
				trb = trb.head(1)
				filtered = filtered.append(trb, ignore_index=True)

		filtered.to_csv(f'{outdir}/filtered.txt', sep='\t')

	elif type == 'BCR':

		data = pd.read_csv(f'{ass_dir}/bracer/filtered_BCR_summary/changeodb.tab', sep='\t')
		data = data[data['FUNCTIONAL'] == True]
		cell_name = set(list(data['CELL']))
		filtered = pd.DataFrame()
		for name in cell_name:
			count_cell = data[data['CELL'] == name]
			count_h = pd.DataFrame(count_cell[count_cell['LOCUS'] == 'H'])
			count_k = pd.DataFrame(count_cell[count_cell['LOCUS'] == 'K'])
			count_l = pd.DataFrame(count_cell[count_cell['LOCUS'] == 'L'])
			count_k_l = count_k.append(count_l)
			if count_h.empty is not True:
				count_h = count_h.sort_values(by='TPM', ascending=False)
				count_h = count_h.head(1)
				filtered = filtered.append(count_h, ignore_index=True)
			if count_k_l.empty is not True:
				count_k_l = count_k_l.sort_values(by='TPM', ascending=False)
				count_k_l = count_k_l.head(1)
				filtered = filtered.append(count_k_l, ignore_index=True)

		filtered.to_csv(f'{outdir}/filtered.txt', sep='\t')

	return filtered


@utils.add_log					
def vdj_sum(args):

	step_name = f"vdj_sum"
	step = Step(args, step_name)

	type = args.type
	ass_dir = args.ass_dir
	sample = args.sample
	outdir = args.outdir
	fastq_dir = args.fastq_dir
	UMI_min = args.UMI_min

	filtered = filtering(type, ass_dir, outdir)

	fqs = glob.glob(f'{fastq_dir}/*.fq')
	matched_bcs = len(fqs)

	stat_file = outdir + '/stat.txt'

	vdj_sum_summary = []
	
	count_umi_file = f'{fastq_dir}/../count.txt'

	count_umi = pd.read_csv(count_umi_file, sep='\t', index_col=0)
	
	all_cells = count_umi.shape[0]

	if type == 'TCR':

		productive_cells = set(filtered['cell_name'].tolist())

		count_umi['mark'] = count_umi['cell_name'].apply(lambda x: "CB" if (x in productive_cells) else "UB")

		count_umi.to_csv(count_umi_file, sep='\t')

		step.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_umi_file))

		productive_cells_num = len(productive_cells)

		TRA_chain = filtered[filtered['locus'] == 'A']
		TRA_chain_num = TRA_chain.shape[0]
		TRB_chain = filtered[filtered['locus'] == 'B']
		TRB_chain_num = TRB_chain.shape[0]

		TRAs, TRBs = [], []
		paired_cell = 0
		for cell in productive_cells:
			tmp1 = TRA_chain[TRA_chain['cell_name'] == cell]
			if tmp1.empty is not True:
				chainA = tmp1['CDR3aa'].tolist()[0]
				TRAs.append(chainA)
			else:
				TRAs.append('NaN')
			
			tmp2 = TRB_chain[TRB_chain['cell_name'] == cell]
			if tmp2.empty is not True:
				chainB = tmp2['CDR3aa'].tolist()[0]
				TRBs.append(chainB)
			else:
				TRBs.append('NaN')
			
			if not tmp1.empty and not tmp2.empty:
				paired_cell += 1

		clonetypes_table = pd.DataFrame()
		clonetypes_table['TRA_chain'] = TRAs
		clonetypes_table['TRB_chain'] = TRBs
		clonetypes_table['Frequency'] = ''

		clonetypes = clonetypes_table.groupby(['TRA_chain', 'TRB_chain']).agg({'Frequency': 'count'})

		sum = clonetypes['Frequency'].sum()
		proportions = []
		for f in list(clonetypes['Frequency']):
			p = f/sum
			p = round(p, 4)
			p = str(p * 100) + '%'
			proportions.append(p)
		clonetypes['Proportion'] = proportions
		clonetypes = clonetypes.sort_values(by='Frequency', ascending=False)
		clonetypes = clonetypes.reset_index()

		clonetypes['clonetypeId'] = [i for i in range(1, (clonetypes.shape[0]+1))]
		clonetypes = clonetypes.reindex(columns=list(['clonetypeId', 'TRA_chain', 'TRB_chain', 'Frequency', 'Proportion']))

		clonetypes.to_csv(f'{outdir}/clonetypes.txt', sep='\t')

		vdj_sum_summary.append({
			'item': 'Estimated Number of Cells',
			'count': productive_cells_num,
			'total_count': all_cells,
		})

		vdj_sum_summary.append({
			'item': 'Cells with TRA',
			'count': TRA_chain_num,
			'total_count': all_cells,
		})

		vdj_sum_summary.append({
			'item': 'Cells with TRB',
			'count': TRB_chain_num,
			'total_count': all_cells,
		})

		vdj_sum_summary.append({
			'item': 'Cells with paired TRA and TRB',
			'count': paired_cell,
			'total_count': all_cells,
		})

		with open(f'{ass_dir}/tmp.txt', 'r') as f:
			medians = []
			for line in f:
				line = line.rstrip('\n').split(':')
				medians.append(int(line[1]))

			vdj_sum_summary.append({
				'item': 'Median UMIs per cell',
				'count': medians[0],
				'total_count': np.nan
			})

			vdj_sum_summary.append({
				'item': 'Median TRA UMIs per cell',
				'count': medians[1],
				'total_count': np.nan	
			})

			vdj_sum_summary.append({
				'item': 'Median TRB UMIs per cell',
				'count': medians[2],
				'total_count': np.nan
			})


	elif type == 'BCR':

		productive_cells = set(filtered['CELL'].tolist())

		productive_cells_num = len(productive_cells)

		count_umi['mark'] = count_umi['cell_name'].apply(lambda x: "CB" if (x in productive_cells) else "UB")

		count_umi.to_csv(count_umi_file, sep='\t')		

		step.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_umi_file))

		filtered_h = filtered[filtered['LOCUS'] == 'H']
		filtered_k = filtered[filtered['LOCUS'] == 'K']
		filtered_l = filtered[filtered['LOCUS'] == 'L']
		filtered_h_count = filtered_h.shape[0]
		filtered_k_count = filtered_k.shape[0]
		filtered_l_count = filtered_l.shape[0]

		IGHs, IGKs, IGLs = [], [], []

		paired_k, paired_l = 0, 0

		for cell in productive_cells:
			tmp1 = filtered_h[filtered_h['CELL'] == cell]
			if tmp1.empty is not True:
				seq = tmp1['JUNCTION'].tolist()[0]
				seq = Seq(seq)
				aaseq = seq.translate()
				IGHs.append(aaseq)
			else:
				IGHs.append('NaN')

			tmp2 = filtered_l[filtered_l['CELL'] == cell]
			if tmp2.empty is not True:
				seq = tmp2['JUNCTION'].tolist()[0]
				seq = Seq(seq)
				aaseq = seq.translate()
				IGLs.append(aaseq)
			else:
				IGLs.append('NaN')

			tmp3 = filtered_k[filtered_k['CELL'] == cell]
			if tmp3.empty is not True:
				seq = tmp3['JUNCTION'].tolist()[0]
				seq = Seq(seq)
				aaseq = seq.translate()
				IGKs.append(aaseq)
			else:
				IGKs.append('NaN')

			if not tmp1.empty and not tmp2.empty:
				paired_l += 1
			if not tmp1.empty and not tmp3.empty:
				paired_k += 1

		clonetypes_table = pd.DataFrame()

		clonetypes_table['IGH_chain'] = IGHs
		clonetypes_table['IGL_chain'] = IGLs
		clonetypes_table['IGK_chain'] = IGKs
		clonetypes_table['Frequency'] = ''

		clonetypes = clonetypes_table.groupby(['IGH_chain', 'IGL_chain', 'IGK_chain']).agg({'Frequency': 'count'})

		Proportion = []
		sum = clonetypes['Frequency'].sum()
		for f in list(clonetypes['Frequency']):
			p = f/sum
			p = round(p, 4)
			p = str(p*100) + '%'
			Proportion.append(p)
		clonetypes['Proportion'] = Proportion
		clonetypes = clonetypes.sort_values(by='Frequency', ascending=False)
		clonetypes = clonetypes.reset_index()

		clonetypes['clonetypeId'] = [i for i in range(1, (clonetypes.shape[0]+1))]
		clonetypes = clonetypes.reindex(columns=list(['clonetypeId', 'IGH_chain', 'IGL_chain', 'IGK_chain', 'Frequency', 'Proportion']))
		clonetypes.to_csv(f'{outdir}/clonetypes.tsv', sep='\t')


		vdj_sum_summary.append({
				'item': 'Estimated Number of Cells',
				'count': productive_cells_num,
				'total_count': all_cells
		})

		vdj_sum_summary.append({
				'item': 'Cells with IGH',
				'count': filtered_h_count,
				'total_count': all_cells
		})	

		vdj_sum_summary.append({
				'item': 'Cells with IGK',
				'count': filtered_k_count,
				'total_count': all_cells
		})

		vdj_sum_summary.append({
				'item': 'Cells with IGL',
				'count': filtered_l_count,
				'total_count': all_cells
		})			

		vdj_sum_summary.append({
				'item': 'Cells with IGH and IGK',
				'count': paired_k,
				'total_count': all_cells
		})

		vdj_sum_summary.append({
				'item': 'Cells with IGH and IGL',
				'count': paired_l,
				'total_count': all_cells
		})

		with open(f'{ass_dir}/tmp.txt', 'r') as f:
			medians=[]
			for line in f:
				line = line.strip('\n').split(':')
				medians.append(int(line[1]))

			vdj_sum_summary.append({
				'item': 'Median UMIs per cell',
				'count': medians[0],
				'total_count': np.nan
			})

			vdj_sum_summary.append({
				'item': 'Median IGH UMIs per cell',
				'count': medians[1],
				'total_count': np.nan
			})

			vdj_sum_summary.append({
				'item': 'Median IGK UMIs per cell',
				'count': medians[2],
				'total_count': np.nan
			})

			vdj_sum_summary.append({
				'item': 'Median IGL UMIs per cell',
				'count': medians[3],
				'total_count': np.nan
			})

	df = pd.DataFrame(vdj_sum_summary, 
		columns=['item', 'count', 'total_count'])

	df['count'] = df['count'].apply(int)
	
	df['percent'] = df['count']/(df.total_count.astype('float')) * 100

	df['percent'] = df['percent'].apply(
		lambda x: round(x, 2)
	)
	df['count'] = df['count'].apply(utils.format_number)


	def percent_str_func(row):
		need_percent = bool(
			re.search("Cells with", row["item"], flags=re.IGNORECASE))
		if need_percent:
			return "(" + str(row["percent"]) + "%)"
		else:
			return ""	

	df['percent_str'] = df.apply(
		lambda row: percent_str_func(row), axis=1
	)	

	def gen_stat(summary, stat_file):
		stat = summary
		stat["new_count"] = stat["count"].astype(str) + stat["percent_str"]
		stat = stat.loc[:, ["item", "new_count"]]
		stat.to_csv(stat_file, sep=":", header=None, index=False)

	gen_stat(df, stat_file)

# clonetype table

	title = 'Clonetypes'
	table_dict = step.get_table(title, 'clonetypes_table', clonetypes)

	step.add_data_item(table_dict=table_dict)

	step.clean_up()


def get_opts_vdj_sum(parser, sub_program):
	if sub_program:
		parser = s_common(parser)
		parser.add_argument('--ass_dir', help='assemble dir', required=True)
		parser.add_argument('--fastq_dir', help='dir contains fastq', required=True)
	parser.add_argument('--type', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
	parser.add_argument('--UMI_min', help='int, min UMI per cell, if not set, will be counted by UMI rank 20', default='auto')



