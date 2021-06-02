import pysam
from collections import defaultdict
import os
import argparse
import datetime
import pandas as pd
from Bio.Seq import Seq
import glob
from celescope.tools import utils
from celescope.tools.Step import Step, s_common
import glob



def tpm_count(ass_dir):
	rec = pd.read_csv(f'{ass_dir}/tracer/filtered_TCRAB_summary/recom' # ass_dir outdir/sample/04.go_assemble
					f'binants.txt', sep='\t')
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


def res_sum(type, ass_dir, outdir):
	filtered = filtering(type, ass_dir, outdir)

	if type == 'TCR':
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
		clone_count = pd.DataFrame(per_count_data['cdr3s_aa'].value_counts())
		clone_count.columns = ["frequency"]
		proportation = []
		sum = clone_count['frequency'].sum()
		for f in list(clone_count['frequency']):
			p = f/sum
			proportation.append(p)
		clone_count['proportation'] = proportation
		clone_count = clone_count.reset_index()
		clone_count.rename(columns={'index': 'cdr3s_aa'}, inplace=True)
		clone_count.to_csv(f'{outdir}/clone_count.tsv', sep='\t')

		return productive_cells, count_a, count_b, paired_cell

	elif type == 'BCR':
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
		clone_count = pd.DataFrame(clones['cdr3s_aa'].value_counts())
		clone_count.columns = ["frequency"]
		proportation = []
		sum = clone_count['frequency'].sum()
		for f in list(clone_count['frequency']):
			p = f/sum
			proportation.append(p)
		clone_count['proportation'] = proportation
		clone_count = clone_count.reset_index()
		clone_count.rename(columns={'index': 'cdr3s_aa'}, inplace=True)
		clone_count.to_csv(f'{outdir}/clone_count.tsv', sep='\t')

		return productive_cells, filtered_h_count, filtered_k_count, filtered_l_count, paired_k, paired_l



def get_stat(fastq_dir, ass_dir, outdir, type):
	fqs = glob.glob(f'{fastq_dir}/*.fq')
	matched_bcs = len(fqs)

	stat_file = outdir + '/stat.txt'
	if type == 'TCR':
		productive_cells, TRA_num, TRB_num, paired_num = res_sum(type, ass_dir, outdir)

		stat_text = pd.DataFrame({
			'item': ['Matched cells', 'Productive cells', 'Cells with TRA', 'Cells with TRB', 'Cells with paired TRA and TRB'], 
			'count': [matched_bcs, productive_cells, TRA_num, TRB_num, paired_num]
		}, 
		columns=['item', 'count'])
		stat_text.to_csv(stat_file, sep=':', header=None, index=False)

	elif type == 'BCR':
		productive_cells, H_num, K_num, L_num, H_K_num, H_L_num = res_sum(type, ass_dir,outdir)

		stat_text = pd.DataFrame({
			'item': ['Matched cells', 'Productive cells', 'Cells with IGH', 'Cells with IGK', 'Cells with IGL', 'Cells with IGH and IGK', 'Cells with IGH and IGL'],
			'count': [matched_bcs, productive_cells, H_num, K_num, L_num, H_K_num, H_L_num]
		}, 
		columns=['item', 'count'])
		stat_text.to_csv(stat_file, sep=":", header=None, index=False)


@utils.add_log					
def vdj_sum(args):
	type = args.type
	ass_dir = args.ass_dir
	sample = args.sample
	outdir = args.outdir
	fastq_dir = args.fastq_dir

	get_stat(fastq_dir, ass_dir, outdir, type)


def get_opts_vdj_sum(parser, sub_program):
	if sub_program:
		parser = s_common(parser)
		parser.add_argument('--ass_dir', help='assemble dir', required=True)
		parser.add_argument('--fastq_dir', help='dir contains fastq', required=True)
	parser.add_argument('--type', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)



