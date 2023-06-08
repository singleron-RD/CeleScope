import pandas as pd
import json
import pysam
import glob
import os
import io
import numpy as np
from collections import defaultdict

from celescope.tools import utils
from celescope.tools.step import s_common, Step
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.flv_CR.match import gen_vj_annotation_metrics
from jinja2 import Environment, FileSystemLoader, select_autoescape


env = Environment(
        loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
        autoescape=select_autoescape(['html', 'xml'])
    )


def get_cltypes(df_contig):
    """
    有一个pandas dataframe,名称为df_contig, 包含如下column：['barcode','chain','cdr3','cdr3_nt']，进行如下操作：
    groupby barcode, 生成一个新的dataframe叫做df_chain，包含新的两列：'cdr3_aa'是{chain}:{cdr3},多个值用分号分隔；‘'cdr3_nt’是{chain}:{cdr3_nt},多个值用分号分隔。
    将df_chain按照cdr3_aa进行groupby, agg生成一个新的dataframe叫做df_cl,barcode这一列进行计数，这一列叫做'frequency'。cdr3_nt这一列如果相同，则保留一个，如果不同，则用分号连接。
    将df_cl按照frequency降序排列。
    """
    grouped = df_contig.groupby('barcode')

    def concatenate_values(x):
        cdr3_aa = ';'.join(f"{chain}:{cdr3}" for chain, cdr3 in zip(x['chain'], x['cdr3']))
        cdr3_nt = ';'.join(f"{chain}:{cdr3_nt}" for chain, cdr3_nt in zip(x['chain'], x['cdr3_nt']))
        return pd.Series({'cdr3s_aa': cdr3_aa, 'cdr3s_nt': cdr3_nt})

    df_chain = grouped.apply(concatenate_values).reset_index()
    barcode_cdr3s_aa_dic = {}
    for _, row in df_chain.iterrows():
        barcode_cdr3s_aa_dic[row['barcode']] = row['cdr3s_aa']

    df_cl = df_chain.groupby('cdr3s_aa').agg({'barcode': 'count', 'cdr3s_nt': lambda x: ';'.join(x.unique())}).reset_index()
    df_cl.rename(columns={'barcode': 'frequency'},inplace=True)
    df_cl = df_cl.sort_values('frequency', ascending=False)

    return df_cl, barcode_cdr3s_aa_dic

def merge_cltypes(cl_origin, cl_refine):
    """
    有两个pandas dataframe：cl_origin和cl_refine,cl_origin包含column：[clonotype_id,frequency,cdr3s_aa,cdr3s_nt], cl_refine包含column:[frequency,cdr3s_aa,cdr3s_nt]
    对cl_refine的每一行做如下处理：如果cdr3s_aa在cl_origin的cdr3s_aa中出现，则将这一行的frequency加在cl_origin第一次出现位置的frequency上，然后从cl_refine中去掉这一行
    将cl_origin和cl_refine合并为df_merge
    为df_merge新生成一列proportion, 值为每一行frequency除以frequncy的sum，保留六位小数。
    将df_merge clonotype为NA的行设为"refine_{num}", num从1开始递增
    """
    cl_origin_copy = cl_origin.copy()
    cl_refine_copy = cl_refine.copy()

    # 遍历cl_refine的每一行
    for index, row in cl_refine_copy.iterrows():
        cdr3s_aa = row['cdr3s_aa']
        mask = cl_origin_copy['cdr3s_aa'] == cdr3s_aa
        
        # 如果cdr3s_aa在cl_origin的cdr3s_aa中出现，则将frequency相加
        if mask.any():
            first_occurrence_index = cl_origin_copy.loc[mask].index[0]
            cl_origin_copy.loc[first_occurrence_index, 'frequency'] += row['frequency']
        
        # 从cl_refine中去掉相应行
            cl_refine_copy.drop(index, inplace=True)

    df_merge = pd.concat([cl_origin_copy, cl_refine_copy])
    df_merge.sort_values(['frequency'], ascending=[False], inplace=True)

    # 计算proportion列
    df_merge['proportion'] = df_merge['frequency'] / df_merge['frequency'].sum()
    df_merge['proportion'] = df_merge['proportion'].round(6)


    # 将clonotype_id为NA的行设为"refine_{num}"
    na_rows = df_merge['clonotype_id'].isna()
    df_merge.loc[na_rows, 'clonotype_id'] = df_merge[na_rows].index.to_series().apply(lambda x: f"refine_{x+1}")

    cdr3s_aa_cl_dic = {}
    for _, row in df_merge.iterrows():
        cdr3s_aa_cl_dic[row['cdr3s_aa']] = row['clonotype_id']

    return df_merge,cdr3s_aa_cl_dic


class Filter_noise:
    """ Filter noise barcode.
    For is_cell=False barcodes, keep the highest UMI chain when its UMI >= coeff * second largest UMI. Then only paired-chain barcodes are keeped.

    First sort in descending order according to barcode, chain, umis, so that the next step does not need to be sorted Perform groupby by barcode and chain.
    If each groupby object has only one data, keep it; if it is more than 1, check whether the highest umi is greater than the second highest umi*1.5, if it is satisfied, keep the highest one, if not, then all not reserved Perform groupby by barcode. 
    If the chain number of each barcode is not equal to 2, it will not be retained.
    When retaining, use a set to record the contig_id, and after all data checks are completed, determine which rows need to be retained according to the contig_id.
    """
    
    def __init__(self, df, coeff, seqtype):

        self.df = df.sort_values("umis", ascending=False)
        self.coeff = coeff
        self.seqtype = seqtype
    
    @utils.add_log
    def __call__(self):
        
        df_sorted = self.df.sort_values(['barcode', 'chain', 'umis'], ascending=[False, False, False])

        # 按barcode和chain进行groupby，保留满足条件的数据
        grouped = df_sorted.groupby(['barcode', 'chain'])
        filtered_data = pd.DataFrame(columns=df_sorted.columns)  # 创建一个空的DataFrame来存储过滤后的数据
        contig_id_set = set()  # 用于记录需要保留的contig_id

        for _group_name, group_data in grouped:
            if len(group_data) == 1:
                contig_id_set.add(group_data['contig_id'].values[0])  
            else:
                sorted_umis = group_data['umis'].values
                highest_umi = sorted_umis[0]
                second_highest_umi = sorted_umis[1]
                if highest_umi > second_highest_umi * self.coeff:
                    contig_id_set.add(group_data.iloc[0]['contig_id'])  # 记录保留的contig_id

        filtered_data = self.df[self.df['contig_id'].isin(contig_id_set)]

        grouped_barcode = filtered_data.groupby('barcode')
        final_filtered_data = pd.DataFrame(columns=filtered_data.columns)
        contig_id_set = set()

        for _group_name, group_data in grouped_barcode:
            if len(group_data['chain'].unique()) == 2:
                if self.seqtype == 'BCR' and 'IGH' not in group_data['chain'].unique():
                    continue
                contig_id_set.update(set(group_data['contig_id']))  # 记录保留的contig_id

        # 根据contig_id确定需要保留的行
        final_filtered_data = filtered_data[filtered_data['contig_id'].isin(contig_id_set)]
        final_filtered_data = final_filtered_data.sort_values(['barcode'])

        return final_filtered_data



class Refine_vdj(Step):
    """
    ## Features
    
    - Refine barcodes where "is_cell=False" and have multi productive chains.

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.assemble_out = args.assemble_out
        self.all_contig_anno = pd.read_csv(glob.glob(f"{self.assemble_out}/all_contig_annotations.csv")[0])
        self.filter_contig_anno = pd.read_csv(glob.glob(f"{self.assemble_out}/filtered_contig_annotations.csv")[0])
        self.all_contig_fasta = glob.glob(f"{self.assemble_out}/all_contig.fasta")[0]
        self.filter_contig_fasta = glob.glob(f"{self.assemble_out}/filtered_contig.fasta")[0]
        self.all_bam = glob.glob(f"{self.assemble_out}/all_contig.bam")[0]
        self.cl_types = f"{self.assemble_out}/clonotypes.csv"
        self.df_cl_origin = pd.read_csv(self.cl_types)
        
        self.seqtype = args.seqtype
        self.coeff = float(args.coeff)
        self.match_cell_barcodes = {}
        
        if args.match_dir != 'None':
            self.match_cell_barcodes, _ = utils.get_barcode_from_match_dir(args.match_dir)

        with open(args.barcode_convert_json, 'r') as f:
            self.tenX_sgr = json.load(f)

        # out files
        self.out_cl_file = f'{self.outdir}/clonotypes.csv'


    def run(self):
        refine_cells, df_merge, df_match = self.refine()
        self.Barcode_rank_plot(refine_cells)
        self.render_html(refine_cells, df_merge, df_match)

    
    @utils.add_log
    def refine(self):
        df_productive = self.all_contig_anno[self.all_contig_anno["productive"] == True]
        df_refine = df_productive[df_productive["is_cell"] == False]
        df_refine = df_refine.groupby("barcode").filter(lambda x: (len(x) > 4))
        df_refine = Filter_noise(df_refine, self.coeff, self.seqtype)()

        # clonotypes
        df_cl_refine,barcode_cdr3s_aa_dic = get_cltypes(df_refine)
        df_cl_merge,cdr3s_aa_cl_dic = merge_cltypes(self.df_cl_origin,df_cl_refine)
        df_cl_merge.to_csv(self.out_cl_file, index=False)
        bc_cl_dic = {}
        for bc in barcode_cdr3s_aa_dic:
            bc_cl_dic[bc] = cdr3s_aa_cl_dic[barcode_cdr3s_aa_dic[bc]]        
        
        # add clonotypes
        df_refine['raw_clonotype_id'] = df_refine['barcode'].replace(bc_cl_dic)

        df_merge = pd.concat([self.filter_contig_anno, df_refine])
        refine_cells = set(df_merge.barcode)
        refine_contigs = set(df_merge.contig_id)

        # generate fasta file
        out_filter_fasta = open(f"{self.outdir}/filtered_contig.fasta", 'w')
        out_match_filter_fasta = open(f"{self.outdir}/matched_contig.fasta", 'w')
        with pysam.FastxFile(self.all_contig_fasta) as f:
            for entry in f:
                name = entry.name
                seq = entry.sequence
                attrs = name.split('_')
                if name in refine_contigs:
                    convert_bc = self.tenX_sgr[attrs[0].split('-')[0]]
                    new_name = convert_bc + '_' + attrs[1] + '_' + attrs[2]
                    out_filter_fasta.write(utils.fasta_line(new_name, seq))
                    if convert_bc in self.match_cell_barcodes:
                        out_match_filter_fasta.write(utils.fasta_line(new_name, seq))

        out_filter_fasta.close()
        out_match_filter_fasta.close()

        # generate annotation file
        df_merge['barcode'] = df_merge['barcode'].apply(lambda x: self.tenX_sgr[x.split('-')[0]])
        df_merge['contig_id'] = df_merge['contig_id'].apply(
            lambda x: self.tenX_sgr[x.split('-')[0]] + '_' + x.split('_')[1] + '_' + x.split('_')[2])
        df_merge['is_cell'] = "TRUE"
        df_merge.to_csv(f"{self.outdir}/filtered_contig_annotations.csv", sep=',', index=False)
        df_match = df_merge[df_merge["barcode"].isin(self.match_cell_barcodes)]
        df_match.to_csv(f"{self.outdir}/matched_contig_annotations.csv", sep=',', index=False)

        return refine_cells, df_merge, df_match

    @utils.add_log
    def Barcode_rank_plot(self, refine_cells):
        dic_umi = defaultdict(set)

        with pysam.AlignmentFile(self.all_bam) as fh:
            for read in fh:
                cb = read.get_tag('CB')
                umi = read.get_tag('UB')
                dic_umi[cb].add(umi)

        df_umi = pd.DataFrame()
        df_umi['barcode'] = list(dic_umi.keys())
        df_umi['UMI'] = [len(dic_umi[i]) for i in dic_umi]
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in refine_cells else 'UB')
        df_umi['barcode'] = df_umi['barcode'].apply(lambda x: self.tenX_sgr[x.split('-')[0]])
        df_umi.to_csv(f"{self.outdir}/count.txt", sep='\t', index=False)


    @utils.add_log
    def render_html(self, refine_cells, df_merge, df_match):
        """
        Generate new refine_vdj html report.
        """
        fh = open(f"{self.outdir}/../.data.json")
        data = json.load(fh)
        fh.close()

        """
        Cells metrics
        """
        cells_reads, total_reads = 0, 0
        with pysam.AlignmentFile(self.all_bam) as f:
            for read in f:
                total_reads += 1
                cb = read.get_tag('CB')
                if cb in refine_cells:
                    cells_reads += 1

        data["cells_summary"]["metric_list"][0]["display"] = str(format(len(refine_cells), ','))
        data["cells_summary"]["metric_list"][1]["display"] = f'{round(cells_reads / total_reads * 100, 2)}%'
        data["cells_summary"]["metric_list"][2]["display"] = str(format(total_reads // len(refine_cells), ','))
        data["cells_summary"]["metric_list"][3]["display"] = str(format(cells_reads // len(refine_cells), ','))
        
        if self.seqtype == "TCR":
            data["cells_summary"]["metric_list"][4]["display"] = np.median(df_merge[df_merge["chain"]=="TRA"].umis)
            data["cells_summary"]["metric_list"][5]["display"] = np.median(df_merge[df_merge["chain"]=="TRB"].umis)
        else:
            data["cells_summary"]["metric_list"][4]["display"] = np.median(df_merge[df_merge["chain"]=="IGH"].umis)
            data["cells_summary"]["metric_list"][5]["display"] = np.median(df_merge[df_merge["chain"]=="IGL"].umis)
            data["cells_summary"]["metric_list"][6]["display"] = np.median(df_merge[df_merge["chain"]=="IGK"].umis)
            
        data["cells_summary"]["chart"] = get_plot_elements.plot_barcode_rank(f"{self.outdir}/count.txt", log_uniform=True)


        """
        Annotation metrics
        """
        metrics_dict = gen_vj_annotation_metrics(df_merge, self.seqtype)
        for k, v in metrics_dict.items():
            for i in data["annotation_summary"]["metric_list"]:
                if i["name"] == k:
                    i["display"] = f'{round(v / len(refine_cells) * 100, 2)}%'

        """
        Match metrics
        """
        if not df_match.empty:
            data["match_summary"]["metric_list"][0]["display"] = str(format(len(set(df_match.barcode)), ','))
            metrics_dict = gen_vj_annotation_metrics(df_match, self.seqtype)
            for k, v in metrics_dict.items():
                for i in data["match_summary"]["metric_list"]:
                    if i["name"] == k:
                        i["display"] = f'{round(v / len(set(df_match.barcode)) * 100, 2)}%'


        """
        Clonotypes table
        """
        df_cl = pd.read_csv(self.out_cl_file)
        title = 'Clonetypes'

        table_dict = self.get_table_dict(
            title=title,
            table_id='clonetypes',
            df_table=df_cl[['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion']]
        )
        data["match_summary"]["table_dict"] = table_dict


        raw_report_html = f"{self.outdir}/../.{self.sample}_raw_report.html"
        report_html = f"{self.outdir}/../{self.sample}_report.html"
        os.system(f"mv {report_html} {raw_report_html}")
        
        template = env.get_template(f'html/{self.assay}/base.html')
        with io.open(report_html, 'w', encoding='utf8') as fh:
            html = template.render(data)
            fh.write(html)


def refine_vdj(args):
    if args.not_refine:
        return
    step_name = 'refine_vdj'
    refine_vdj_obj = Refine_vdj(args, step_name)
    refine_vdj_obj.run()


def get_opts_refine_vdj(parser, sub_program):
    parser.add_argument("--not_refine", help="Do not perform the refine step. ", action="store_true")
    parser.add_argument("--coeff",
        help="coefficient will affect auto and snr noise filter, recommend 1.5 for auto", 
        type=float, default=1.5)
    parser.add_argument("--seqtype", help="TCR or BCR", choices=["TCR", "BCR"], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument("--assemble_out", help="directory of cellranger assemble result", required=True)
        parser.add_argument("--match_dir", help="scRNA-seq match directory", required=True)
        parser.add_argument("--barcode_convert_json", help="json file", required=True)
        
    return parser