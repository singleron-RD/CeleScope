import os
import pandas as pd
from Bio.Seq import Seq
import numpy as np
from celescope.tools import utils
from celescope.tools.Step import Step, s_common
from celescope.tools.cellranger3 import get_plot_elements
from celescope.tracer_vdj.go_assemble import percent_str_func, gen_stat


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


def filtering(Seqtype, ass_dir, outdir):

	if not os.path.exists(outdir):
		os.makedirs(outdir)

	if Seqtype == 'TCR':
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

	elif Seqtype == 'BCR':

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

		results = filtering(Seqtype, ass_dir, outdir)

		stat_file = outdir + '/stat.txt'

		vdj_sum_summary = []
		
		count_umi_file = f'{fastq_dir}/../{self.sample}_count.txt'

		count_umi = pd.read_csv(count_umi_file, sep='\t', index_col=0)
		
		all_cells = count_umi.shape[0]

		if Seqtype == 'TCR':

			productive_cells = set(results['cell_name'].tolist())

			count_umi['mark'] = count_umi['cell_name'].apply(lambda x: "CB" if (x in productive_cells) else "UB")

			count_umi.to_csv(count_umi_file, sep='\t')

			self.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_umi_file))

			productive_cells_num = len(productive_cells)

			TRA_chain = results[results['locus'] == 'A']
			TRA_chain_num = TRA_chain.shape[0]
			TRB_chain = results[results['locus'] == 'B']
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

			clonetypes['clonetypeId'] = [i for i in range(1, (clonetypes.shape[0]+1))]
			clonetypes = clonetypes.reindex(columns=list(['clonetypeId', 'TRA_chain', 'TRB_chain', 'Frequency', 'Proportion']))

			clonetypes.to_csv(f'{outdir}/clonetypes.tsv', sep='\t')

			vdj_sum_summary.append({
				'item': 'Estimated Number of Cells',
				'count': productive_cells_num,
				'total_count': all_cells,
			})

			vdj_sum_summary.append({
				'item': 'Cells with TRA',
				'count': TRA_chain_num,
				'total_count': productive_cells_num,
			})

			vdj_sum_summary.append({
				'item': 'Cells with TRB',
				'count': TRB_chain_num,
				'total_count': productive_cells_num,
			})

			vdj_sum_summary.append({
				'item': 'Cells with paired TRA and TRB',
				'count': paired_cell,
				'total_count': productive_cells_num,
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


		elif Seqtype == 'BCR':

			productive_cells = set(results['CELL'].tolist())

			productive_cells_num = len(productive_cells)

			count_umi['mark'] = count_umi['cell_name'].apply(lambda x: "CB" if (x in productive_cells) else "UB")

			count_umi.to_csv(count_umi_file, sep='\t')		

			self.add_data_item(chart=get_plot_elements.plot_barcode_rank(count_umi_file))

			results_h = results[results['LOCUS'] == 'H']
			results_k = results[results['LOCUS'] == 'K']
			results_l = results[results['LOCUS'] == 'L']
			results_h_count = results_h.shape[0]
			results_k_count = results_k.shape[0]
			results_l_count = results_l.shape[0]

			IGHs, IGKs, IGLs = [], [], []

			paired_k, paired_l = 0, 0

			for cell in productive_cells:
				tmp1 = results_h[results_h['CELL'] == cell]
				if tmp1.empty is not True:
					seq = tmp1['JUNCTION'].tolist()[0]
					seq = Seq(seq)
					aaseq = seq.translate()
					IGHs.append(aaseq)
				else:
					IGHs.append('NaN')

				tmp2 = results_l[results_l['CELL'] == cell]
				if tmp2.empty is not True:
					seq = tmp2['JUNCTION'].tolist()[0]
					seq = Seq(seq)
					aaseq = seq.translate()
					IGLs.append(aaseq)
				else:
					IGLs.append('NaN')

				tmp3 = results_k[results_k['CELL'] == cell]
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
					'count': results_h_count,
					'total_count': productive_cells_num
			})	

			vdj_sum_summary.append({
					'item': 'Cells with IGK',
					'count': results_k_count,
					'total_count': productive_cells_num
			})

			vdj_sum_summary.append({
					'item': 'Cells with IGL',
					'count': results_l_count,
					'total_count': productive_cells_num
			})			

			vdj_sum_summary.append({
					'item': 'Cells with paired IGH and IGK',
					'count': paired_k,
					'total_count': productive_cells_num
			})

			vdj_sum_summary.append({
					'item': 'Cells with paired IGH and IGL',
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

		df['percent_str'] = df.apply(
			lambda row: percent_str_func(row), axis=1
		)	

		gen_stat(df, stat_file)

	# clonetype table

		title = 'Clonetypes'
		table_dict = self.get_table(title, 'clonetypes_table', clonetypes)

		self.add_data_item(table_dict=table_dict)

		self.clean_up()

		os.remove(f'{ass_dir}/tmp.txt')


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



