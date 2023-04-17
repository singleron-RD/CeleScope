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
from celescope.tools.plotly_plot import Bar_plot
from celescope.flv_CR.match import gen_clonotypes_table, gen_vj_annotation_metrics
from jinja2 import Environment, FileSystemLoader, select_autoescape


env = Environment(
        loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
        autoescape=select_autoescape(['html', 'xml'])
    )


class Filter_noise:
    """ Filter noise barcode.
    auto:
    keep paired-chain barcode when highest umi of contig >= coeff * second highest umi of contig.
    
    snr:
    https://www.nature.com/articles/s41598-017-18303-z#:~:text=The%20signal%2Dto%2Dnoise%20ratio,of%20background%20pixel%20intensity%2C%20respectively.
    keep paired-chain barcode when SNR value >= coefficient.
    The signal-to-noise ratio (SNR) of the single molecules was defined as SNR = (S − B)/σ,
    where S is the peak single molecule pixel intensity and B and 
    σ are the average and standard deviation of background pixel intensity, respectively.
    
    not_filter:
    for paired-chain barcode, keep highest tra and trb chain.
    """
    
    def __init__(self, df, filter_method, coeff, seqtype):

        self.df = df.sort_values("umis", ascending=False)
        self.method = filter_method
        self.coeff = coeff
        self.seqtype = seqtype
    
    @utils.add_log
    def __call__(self):
        
        if self.method == "not_filter":
            self.df = self.not_filter(self.df, self.seqtype)
        
        else:
            bc_chain_dict = self.df.groupby("barcode")["chain"].apply(set).to_dict()
            if self.seqtype == "TCR":
                pair_chain_dict = {key: value for key, value in bc_chain_dict.items() if len(value)==2}
            else:
                pair_chain_dict = {key: value for key, value in bc_chain_dict.items() if len(value)>=2 and "IGH" in value}
            self.df = self.df[self.df["barcode"].isin(pair_chain_dict)]
            self.df = self.df.sort_values(["barcode","umis"], ascending=[False, False])
            
            # multi heavy or light chains
            df_multi_chain = self.df.groupby("barcode").filter(lambda x: (len(x) > 2))
            # single heavy-light pair chain
            df_pair_chain = self.df.groupby("barcode").filter(lambda x: (len(x) == 2))
            
            if self.seqtype == "TCR":
                df_heavy = df_multi_chain[df_multi_chain['chain'] == 'TRB']
                df_light = df_multi_chain[df_multi_chain['chain'] == 'TRA']
            else:
                df_heavy = df_multi_chain[(df_multi_chain['chain'] == 'IGH')] 
                df_light = df_multi_chain[(df_multi_chain['chain'] == 'IGL') | (df_multi_chain['chain'] =='IGK')]
            light_dict = df_light.groupby("barcode")["umis"].apply(lambda x: x.tolist()).to_dict()
            heavy_dict = df_heavy.groupby("barcode")["umis"].apply(lambda x: x.tolist()).to_dict()

            if self.method == "auto":
                for chain_dict in [light_dict, heavy_dict]:
                    for k in chain_dict:
                        if len(chain_dict[k]) == 1:
                            chain_dict[k].append(0.1)

                # highest umi >= coeff * second highest umi
                filter_noise_light = {key for key,value in light_dict.items() if value[0]/value[1] >= self.coeff}
                filter_noise_heavy = {key for key,value in heavy_dict.items() if value[0]/value[1] >= self.coeff}
                
            elif self.method == "snr":
                for chain_dict in [light_dict, heavy_dict]:
                    for k in chain_dict:
                        if len(chain_dict[k]) == 1:
                            chain_dict[k].append(chain_dict[k][0])
                
                # snr value >= coeff
                filter_noise_light = self.snr_filter(light_dict, self.coeff)
                filter_noise_heavy = self.snr_filter(heavy_dict, self.coeff)
                        
            filter_noise_barcode = filter_noise_light & filter_noise_heavy | set(df_pair_chain.barcode)
            self.df = self.df[self.df["barcode"].isin(filter_noise_barcode)]
        
        return self.df

    @staticmethod
    def not_filter(df, seqtype):
        """do not filter any noise, keep all paired-chain barcodes"""
        
        if seqtype == "TCR":
            df = df.groupby(["barcode", "chain"], as_index=False).head(1)
            df = df.groupby("barcode").filter(lambda x: (len(x) > 1))
        else:
            df_chain_heavy = df[(df['chain'] == 'IGH')]
            df_chain_light = df[(df['chain'] == 'IGL') | (df['chain'] =='IGK')]
            df_chain_heavy.drop_duplicates(['barcode'], inplace=True)
            df_chain_light.drop_duplicates(['barcode'], inplace=True)
            df = pd.concat([df_chain_heavy, df_chain_light])
        
        return df
    
    @staticmethod
    def snr_filter(chain_dict, coeff):
        """ calculate SNR for each barcode
        SNR = (S − B)/σ
        
        :param chain_dict: {'AAACATCGAAGGACACCCTAATCC': [7, 3, 1],
                         'AAACATCGAATCCGTCCATCAAGT': [10, 2, 1],
                         'AAACATCGAATCCGTCCCGACAAC': [4, 3],
                         }
        :return: {AAACATCGAAGGACACCCTAATCC, AAACATCGAATCCGTCCATCAAGT}
        """
        filter_noise_barcode = set()
        for k, v in chain_dict.items():
            S, noise = v[0], v[1:]
            B = np.mean(noise)
            O = np.std(noise)
            if O == 0:
                filter_noise_barcode.add(k)
                continue
            SNR = (S - B) / O
            if SNR >= coeff:
                filter_noise_barcode.add(k)

        return filter_noise_barcode


class Refine_vdj(Step):
    """
    ## Features
    
    - Refine barcodes where "is_cell=False" and have multi productive chains.
    
    - There are three methods to filter noise: SNR, AUTO, NOT_FILTER

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.assemble_out = args.assemble_out
        self.all_contig_anno = pd.read_csv(glob.glob(f"{self.assemble_out}/all_contig_annotations.csv")[0])
        self.filter_contig_anno = pd.read_csv(glob.glob(f"{self.assemble_out}/filtered_contig_annotations.csv")[0])
        self.all_contig_fasta = glob.glob(f"{self.assemble_out}/all_contig.fasta")[0]
        self.filter_contig_fasta = glob.glob(f"{self.assemble_out}/filtered_contig.fasta")[0]
        self.all_bam = glob.glob(f"{self.assemble_out}/all_contig.bam")[0]
        
        self.seqtype = args.seqtype
        self.filter_method = args.filter_method
        self.coeff = float(args.coeff)
        self.match_cell_barcodes = {}
        
        if args.match_dir != 'None':
            self.match_cell_barcodes, _ = utils.get_barcode_from_match_dir(args.match_dir)

        with open(args.barcode_convert_json, 'r') as f:
            self.tenX_sgr = json.load(f)


    def run(self):
        refine_cells, df_merge, df_match = self.refine()
        self.Barcode_rank_plot(refine_cells)
        self.render_html(refine_cells, df_merge, df_match)

    
    @utils.add_log
    def refine(self):
        df_productive = self.all_contig_anno[self.all_contig_anno["productive"] == True]
        df_refine = df_productive[df_productive["is_cell"] == False]
        df_refine = Filter_noise(df_refine, self.filter_method, self.coeff, self.seqtype)()
        
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

            gen_clonotypes_table(df_match, f"{self.outdir}/matched_clonotypes.csv")

        """
        Clonotypes table
        """
        gen_clonotypes_table(df_merge, f"{self.outdir}/clonotypes.csv")
        title = 'Clonetypes'
        raw_clonotypes = pd.read_csv(f"{self.outdir}/clonotypes.csv", sep=',', index_col=None)
        raw_clonotypes['ClonotypeID'] = raw_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        raw_clonotypes['Frequency'] = raw_clonotypes['frequency']
        raw_clonotypes['Proportion'] = raw_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        raw_clonotypes['CDR3_aa'] = raw_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))

        table_dict = self.get_table_dict(
            title=title,
            table_id='clonetypes',
            df_table=raw_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']]
        )
        data["match_summary"]["table_dict"] = table_dict

        raw_clonotypes['ClonotypeID'] = raw_clonotypes['ClonotypeID'].astype("int")
        raw_clonotypes.sort_values(by=['ClonotypeID'], inplace=True)
        Barplot = Bar_plot(df_bar=raw_clonotypes).get_plotly_div()
        data["match_summary"]["Barplot"] = Barplot

        template = env.get_template(f'html/{self.assay}/base.html')
        report_html = f"{self.outdir}/{self.sample}_refine_report.html"
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
    parser.add_argument("--not_refine", help="Run the VDJ refine step. ", action="store_true")
    parser.add_argument("--filter_method", help="filter noise method", choices=["auto", "snr", "not_filter"], default="auto")
    parser.add_argument("--coeff",
                        help="coefficient will affect auto and snr noise filter, recommend 1.5 for auto, 10 for snr",
                        default=1.5
                        )
    parser.add_argument("--seqtype", help="TCR or BCR", choices=["TCR", "BCR"], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument("--assemble_out", help="directory of cellranger assemble result", required=True)
        parser.add_argument("--match_dir", help="scRNA-seq match directory", required=True)
        parser.add_argument("--barcode_convert_json", help="json file", required=True)
        
    return parser