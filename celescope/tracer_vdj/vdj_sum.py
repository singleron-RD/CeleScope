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
	
	count_umi = f'{fastq_dir}/../umi_count.tsv'

	if type == 'TCR':

		step.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_umi))

		count_a = filtered[filtered['locus'] == 'A'].shape[0]
		count_b = filtered[filtered['locus'] == 'B'].shape[0]
		paired_cell = pd.DataFrame(filtered['cell_name'].value_counts())
		productive_cells = paired_cell.shape[0]
		unpaired_cell = paired_cell[paired_cell['cell_name'] == 1]
		paired_cell = paired_cell[paired_cell['cell_name'] == 2]
		paired_cell = list(paired_cell.index)

		aaseqs = []
		for cell in paired_cell:
			temp = filtered[filtered['cell_name'] == cell]
			temp_loci = list(temp['locus'])
			temp_aaseq = list(temp['CDR3aa'])
			string = 'TR{}:C{}F;TR{}:C{}F'.format(temp_loci[0], temp_aaseq[0], temp_loci[1], temp_aaseq[1])
			aaseqs.append(string)

		for cell in list(unpaired_cell.index):
			temp = filtered[filtered['cell_name'] == cell]
			temp_loci = list(temp['locus'])
			temp_aaseq = list(temp['CDR3aa'])
			string = 'TR{}:C{}F'.format(temp_loci[0], temp_aaseq[0])
			aaseqs.append(string)

		per_count_data = pd.DataFrame()
		per_count_data['cdr3s_aa'] = aaseqs
		clonetypes = pd.DataFrame(per_count_data['cdr3s_aa'].value_counts())
		clonetypes.columns = ["Frequency"]
		Percent = []
		sum = clonetypes['Frequency'].sum()
		for f in list(clonetypes['Frequency']):
			p = f/sum
			Percent.append(p)
		clonetypes['Percent'] = Percent
		clonetypes = clonetypes.reset_index()
		clonetypes.rename(columns={'index': 'cdr3s_aa'}, inplace=True)
		clonetypes.to_csv(f'{outdir}/clonetypes.tsv', sep='\t')

		vdj_sum_summary.append({
			'item': 'Estimated Number of Cells',
			'count': matched_bcs,
			'total_count': matched_bcs,
		})

		vdj_sum_summary.append({
			'item': 'Productive cells',
			'count': productive_cells,
			'total_count': matched_bcs
		})

		vdj_sum_summary.append({
			'item': 'Cells with TRA',
			'count': count_a,
			'total_count': matched_bcs,
		})

		vdj_sum_summary.append({
			'item': 'Cells with TRB',
			'count': count_b,
			'total_count': matched_bcs,
		})

		vdj_sum_summary.append({
			'item': 'Cells with paired TRA and TRB',
			'count': len(paired_cell),
			'total_count': matched_bcs,
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

		step.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_umi))

		filtered_h = filtered[filtered['LOCUS'] == 'H']
		filtered_k = filtered[filtered['LOCUS'] == 'K']
		filtered_l = filtered[filtered['LOCUS'] == 'L']
		filtered_h_count = filtered_h.shape[0]
		filtered_k_count = filtered_k.shape[0]
		filtered_l_count = filtered_l.shape[0]

		paired_cell = pd.DataFrame(filtered['CELL'].value_counts())
		productive_cells = paired_cell.shape[0]	

		paired_cell = pd.DataFrame(filtered['CELL'].value_counts())
		productive_cells = paired_cell.shape[0]
		unpaired_cell = paired_cell[paired_cell['CELL'] == 1]
		paired_cell = paired_cell[paired_cell['CELL'] == 2]
		paired_k = 0
		paired_l = 0

		clones = pd.DataFrame()
		cells = list(paired_cell.index)
		aaseqs = []

		for cell in cells:
			if 'K' in list(filtered[filtered['CELL'] == cell]['LOCUS']):
				paired_k += 1
			elif 'L' in list(filtered[filtered['CELL'] == cell]['LOCUS']):
				paired_l += 1
			tep = filtered[filtered['CELL'] == cell]
			tep_loci = list(tep['LOCUS'])
			cdr3 = list(tep['JUNCTION'])
			aaseq = []
			for seq in cdr3:
				seq = Seq(seq)
				seq = seq.translate()
				aaseq.append(seq)
			string = 'IG{}:{};IG{}:{}'.format(tep_loci[0], aaseq[0], tep_loci[1], aaseq[1])
			aaseqs.append(string)

		for cell in list(unpaired_cell.index):
			cells.append(cell)
			locus = list(filtered[filtered['CELL'] == cell]['LOCUS'])
			cdr3 = list(filtered[filtered['CELL'] == cell]['JUNCTION'])
			seq = Seq(cdr3[0])
			seq = seq.translate()
			string = 'IG{}:{}'.format(locus[0], seq)
			aaseqs.append(string)

		clones['CELLS'] = cells

		clones["cdr3s_aa"] = aaseqs
		clonetypes = pd.DataFrame(clones['cdr3s_aa'].value_counts())
		clonetypes.columns = ["Frequency"]
		Percent = []
		sum = clonetypes['Frequency'].sum()
		for f in list(clonetypes['Frequency']):
			p = f/sum
			Percent.append(p)
		clonetypes['Percent'] = Percent
		clonetypes = clonetypes.reset_index()
		clonetypes.rename(columns={'index': 'cdr3s_aa'}, inplace=True)
		clonetypes.to_csv(f'{outdir}/clonetypes.tsv', sep='\t')


		vdj_sum_summary.append({
				'item': 'Matched cells',
				'count': matched_bcs,
				'total_count': matched_bcs
		})

		vdj_sum_summary.append({
				'item': 'Productive cells',
				'count': productive_cells,
				'total_count': matched_bcs
		})

		vdj_sum_summary.append({
				'item': 'Cells with IGH',
				'count': filtered_h_count,
				'total_count': matched_bcs
		})	

		vdj_sum_summary.append({
				'item': 'Cells with IGK',
				'count': filtered_k_count,
				'total_count': matched_bcs
		})

		vdj_sum_summary.append({
				'item': 'Cells with IGL',
				'count': filtered_l_count,
				'total_count': matched_bcs
		})			

		vdj_sum_summary.append({
				'item': 'Cells with IGH and IGK',
				'count': paired_k,
				'total_count': matched_bcs
		})

		vdj_sum_summary.append({
				'item': 'Cells with IGH and IGL',
				'count': paired_l,
				'total_count': matched_bcs
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

	clonetypes['Percent'] = clonetypes['Percent'].apply(lambda x: str(round(x*100, 2)) + '%')
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



