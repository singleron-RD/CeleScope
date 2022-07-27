import pandas as pd
import pysam
import copy
import subprocess
import json

from celescope.tools import analysis_wrapper
from collections import defaultdict
from celescope.tools import utils
from celescope.tools.capture.threshold import Auto
from celescope.tools.step import Step, s_common
from celescope.flv_trust4.__init__ import CHAIN, PAIRED_CHAIN, TOOLS_DIR
from celescope.tools.emptydrop_cr import get_plot_elements


ASSEMBLE_SUFFIX = 'assembled_reads.fa'
ANNOTATE_SUFFIX = 'annotate.fa'
TRUST_REPORT_SUFFIX = 'filter_report.tsv'
BARCODE_REPORT_SUFFIX = 'barcode_report.tsv'
BARCODE_REPORT_FILTER_SUFFIX = 'barcode_filter_report.tsv'


def target_cell_calling(df_UMI_sum, expected_target_cell_num=3000, target_barcodes=None, weight=6, coef=5, 
    percentile=85, umi_col='umis'):
    """
    Args:
        df_UMI_sum: A dataframe with columns highest umi's contig and UMI.
    
    Returns:
        target_contigs_id: list
    >>> df_UMI_sum = pd.DataFrame({"contig_id": ["A", "B", "C", "D", "E"], "UMI": [1, 2, 1, 30, 40]})
    >>> target_contigs_id = target_cell_calling(df_UMI_sum, expected_target_cell_num=5, percentile=80, coef=5, target_barcodes=["A", "C"])
    >>> target_contigs_id == {'A_1', 'C_1', 'D_1', 'E_1'}
    True
    """
    if target_barcodes != None:
        target_barcodes = {i for i in target_barcodes}
    umi_threshold = Auto(list(df_UMI_sum[umi_col]), expected_cell_num=expected_target_cell_num, coef=coef, percentile=percentile).run()

    # avoid change the original dataframe
    df_temp = df_UMI_sum.copy()
    if target_barcodes:
        df_temp[umi_col] = df_temp.apply(
            lambda row:  row[umi_col] * weight if row['barcode'] in target_barcodes else row[umi_col], axis=1)
             
    target_contigs = set(df_temp.loc[df_temp[umi_col] >= umi_threshold].contig_id)

    return target_contigs


class Summarize(Step):
    """
    ## Features
    - CDR3 filtering: contain stop condon, length <=5, etc..

    - If barcode A's two chains CDR3s are identical to another barcode B, and A's chain abundance is significantly lower than B's, filter A.

    - If `--target_cell_barcode` is provided, the UMI counts of all contigs originated from target cells are multiplied by a weight(default: 6.0) to better distinguish signal from background noise. `--target_cell_barcode` comes from the cell type annotation results of the RNA library.

    - Cell-calling is similar to the rna cell-calling algorithm.

    ## Output
    - `04.summarize/clonetypes.tsv` High-level description for each clonotype.

    - `04.summarize/{sample}_all_contig.csv` High-level and detailed annotation for each contig.

    - `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.

    - `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of {sample}_all_contig.csv.
    
    - `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.ref = args.ref
        self.fq2 = args.fq2
        self.diffuseFrac = args.diffuseFrac
        self.assembled_fa = f'{args.assemble_out}/{self.sample}_{ASSEMBLE_SUFFIX}'
        self.trust_report = f'{args.assemble_out}/{self.sample}_{TRUST_REPORT_SUFFIX}'
        self.annot = f'{args.assemble_out}/{self.sample}_{ANNOTATE_SUFFIX}'

        self.match_dict = utils.parse_match_dir(args.match_dir)
        self.chains, self.paired_groups = self._parse_seqtype(self.seqtype)

        # if --diffuseFrac provided
        if self.diffuseFrac:
            self.barcode_report = f'{args.assemble_out}/{self.sample}_{BARCODE_REPORT_FILTER_SUFFIX}'
        else:
            self.barcode_report = f'{args.assemble_out}/{self.sample}_{BARCODE_REPORT_SUFFIX}'
        
        self.coef = int(args.coef)
        self.target_weight = args.target_weight

        """
        There are three method. 
        1. --expected_target_cell_num, Expected assembled T or B cell number. eg 3000
        2. --target_cell_barcode Auto, Based on auto assigned cell type annotation results of the RNA library.
        3. --target_cell_barcode Absolute path of a plain text file with one barcode per line.(Recommend, need manual-assign first)
        """
        if not args.target_cell_barcode:
            self.target_barcodes = None
            self.expected_target_cell_num = args.expected_target_cell_num
        elif args.target_cell_barcode == 'Auto':
            self.target_barcodes, self.expected_target_cell_num = self._auto_assign()
        else:
            self.target_barcodes, self.expected_target_cell_num = utils.read_one_col(args.target_cell_barcode)
    
    @staticmethod
    def _parse_seqtype(seqtype):
        """Parse BCR or TCR

        :param seqtype
        :return CHAIN[seqtype]: 'TCR': ['TRA', 'TRB'], 'BCR': ['IGH', 'IGL', 'IGK']
        	    PAIRED_CHAIN[seqtype]: ['TRA_TRB'], 'BCR': ['IGH_IGL', 'IGH_IGK']
        """
        return CHAIN[seqtype], PAIRED_CHAIN[seqtype]
    
    @staticmethod
    @utils.add_log
    def convert_barcode_report(barcode_report, outdir):
        cmd = (
            f'{TOOLS_DIR}/trust-barcoderep-to-10X.pl '
            f'{barcode_report} '
            f'{outdir} '
        )
        Summarize.convert_barcode_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
    
    def add_cell_num_metric(self, df, name):
        """
        add cell number after each filtering step to .metrics.json
        do not show in HTML report
        """
        cell_num = len(set(df[df['productive']==True].barcode))
        self.add_metric(name, cell_num, show=False)
        return cell_num
    
    @utils.add_log
    def _auto_assign(self):
        """
        --target_cell_barcode Auto
        Auto assign cell type of the RNA library based on marker gene file.
        target_cell_barcode requires >=2 identified marker genes of T/B cells in top10 marker genes for each cluster.
        Return target_barcodes, expected_target_cell_num.
        """
        species, cell_type = Summarize.get_cell_species(self.ref, self.seqtype)

        with open(f'{TOOLS_DIR}/Immune_marker.json', 'r') as f:
            marker_dict = json.load(f)

        report_runner = analysis_wrapper.Report_runner(self.args)
        df_tsne, df_marker = report_runner.get_df()
        
        bc_cluster_dict = df_tsne.reset_index().groupby('cluster')['barcode'].apply(lambda x: x.tolist()).to_dict()
        df_marker_top5 = df_marker.groupby('cluster').head(5)
        df_marker_top5['celltype'] = None

        target_genes = set(df_marker_top5.gene)
        for gene in target_genes:
            if gene in marker_dict[species][cell_type]:
                df_marker_top5.loc[df_marker_top5['gene']==gene, 'celltype'] = cell_type

        target_clusters = df_marker_top5[df_marker_top5['celltype']==cell_type].cluster.tolist()
        target_clusters = [int(x.split(' ')[-1]) for x in target_clusters]
        target_clusters = set([x for x in target_clusters if target_clusters.count(x) > 1])
        
        if not target_clusters:
            return None, self.args.expected_target_cell_num

        target_barcodes = []
        for cluster in target_clusters:
            target_barcodes += bc_cluster_dict[cluster]

        return target_barcodes, len(target_barcodes)

    @staticmethod
    def get_cell_species(ref, seqtype):

        if ref == 'GRCm38':
            species = 'mouse'
        else:
            species = 'human'

        if seqtype =='BCR':
            target_cell_type = 'B_cells'
        else:
            target_cell_type = 'T_cells'

        return species, target_cell_type

    @utils.add_log
    def parse_contig_file(self):
        """
        Convert to all_contig_annotation file in 10X format.
        Generate all_contig_fasta file.
        """
        self.convert_barcode_report(self.barcode_report, outdir=f'{self.outdir}/{self.sample}')

        if self.seqtype == 'BCR':
            df = pd.read_csv(f'{self.outdir}/{self.sample}_b.csv')
        else:
            df = pd.read_csv(f'{self.outdir}/{self.sample}_t.csv')

        df['productive'] = df['full_length']
        contig_set = set(df.contig_id)

        # generate all contig fasta file
        # add length of each contig. 
        len_dict = dict()
        all_fa = open(f'{self.outdir}/{self.sample}_all_contig.fasta','w')
        with pysam.FastxFile(self.annot) as fa:
            for read in fa:
                len_dict[read.name] = read.comment.split(' ')[0]
                if read.name in contig_set:
                    sequence = read.sequence
                    all_fa.write('>' + read.name + '\n' + sequence + '\n')    
        all_fa.close()
        df['length'] = df['contig_id'].apply(len_dict.get)
        
        return df

    @utils.add_log
    def cell_calling(self, df):
        """
        Common filtering based on CDR3:
        Filter nonfunctional CDR3(shown 'out_of_frame' in cdr3 report), or CDR3 sequences containing "N" in the nucleotide sequence.
        Keep CDR3aa start with C.
        Keep CDR3aa length >= 5.
        Keep no stop codon in CDR3aa.
        Filter low abundance contigs based on a umi cut-off. 

        DiffuseFrac filtering(option --diffuseFrac needed):
        If cell A's two chains CDR3s are identical to another cell B, 
        and A's chain abundance is significantly lower than B's (--diffuseFrac), filter A.

        Target cell barcodes filtering(option, --target_cell_barcode needed):
        Filter low abundance contigs based on a umi cut-off.
        The umi counts of all contigs originated from B cells are multiplied by a weight to 
        better distinguish signal from background noise.
        """
        df.sort_values(by='umis', ascending=False, inplace=True)
        if self.seqtype == 'BCR':
            df_chain_heavy = df[df['chain']=='IGH']
            df_chain_light = df[(df['chain']=='IGK') | (df['chain']=='IGL')]
        else:
            df_chain_heavy = df[df['chain'] == 'TRA']
            df_chain_light = df[df['chain'] == 'TRB']
        df_chain_heavy = df_chain_heavy.drop_duplicates(['barcode'])
        df_chain_light = df_chain_light.drop_duplicates(['barcode'])
        df_for_clono = pd.concat([df_chain_heavy, df_chain_light], ignore_index=True)
        
        self.add_cell_num_metric(df_for_clono, 'Total Assembled Cell Number')
        
        # Common filtering
        trust_report = pd.read_csv(self.trust_report, sep='\t')
        correct_cdr3 = set(df_for_clono.cdr3).intersection(set(trust_report.CDR3aa))
        correct_cdr3 = [i for i in correct_cdr3 if i.startswith('C')]
        correct_cdr3 = [i for i in correct_cdr3 if len(i)>=5]
        correct_cdr3 = [i for i in correct_cdr3 if 'UAG' or 'UAA' or 'UGA' not in i]
        df_for_clono = df_for_clono[df_for_clono['cdr3'].isin(correct_cdr3)]

        self.add_cell_num_metric(df_for_clono, 'Cell Number after CDR3 filtering')
        
        # Filter low abundance contigs based on a umi cut-off
        if self.seqtype == 'BCR':
            df_chain_heavy = df_for_clono[df_for_clono['chain']=='IGH']
            df_chain_light = df_for_clono[(df_for_clono['chain']=='IGK') | (df_for_clono['chain']=='IGL')]
        else:
            df_chain_heavy = df_for_clono[df_for_clono['chain'] == 'TRA']
            df_chain_light = df_for_clono[df_for_clono['chain'] == 'TRB']

        filtered_congtigs_id = set()
        for _df in [df_chain_heavy, df_chain_light]:
            target_contigs = target_cell_calling(
            _df, 
            expected_target_cell_num=self.expected_target_cell_num, 
            target_barcodes=self.target_barcodes,
            weight = self.target_weight,
            coef = self.coef
            )
            filtered_congtigs_id = filtered_congtigs_id | target_contigs       
        
        df_for_clono = df_for_clono[df_for_clono.contig_id.isin(filtered_congtigs_id)]
        self.add_cell_num_metric(df_for_clono, 'Cell Number after UMI filtering')
        
        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        cell_barcodes, filtered_contig = set(df_for_clono_pro['barcode']), set(df_for_clono_pro['contig_id'])


        return df_for_clono, cell_barcodes, filtered_contig

    @utils.add_log
    def filter_fasta(self, cell_barcodes):
        """Filter all contig fasta file by barcodes which are identified to be cell.

        :param cell_barcodes: all barcodes identified to be cell.
        """
        all_contig_fasta = f'{self.outdir}/{self.sample}_all_contig.fasta'
        filter_contig_fasta = f'{self.outdir}/{self.sample}_filtered_contig.fasta'

        filter_contig_fasta = open(filter_contig_fasta,'w')
        with pysam.FastxFile(all_contig_fasta) as fa:
            for read in fa:
                name = read.name
                cb = name.split('_')[0]
                sequence = read.sequence
                if cb in cell_barcodes:
                    filter_contig_fasta.write('>' + name + '\n' + sequence + '\n')

        filter_contig_fasta.close()

    @utils.add_log
    def parse_clonotypes(self, df, df_for_clono, cell_barcodes, filtered_contig):
        """Parse clonotypes from CDR3 and manually add clonotype id for each contig.

        :param df: original contig file.
        :param df_for_clono: contig info after filter.
        :param cell_barcodes: all barcodes identified to be cell.
        :return df_filter_contig: filtered contigs by cell barcodes.
        """
        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True].copy()
        df_for_clono_pro['chain_cdr3aa'] = df_for_clono_pro.loc[:, ['chain', 'cdr3']].apply(':'.join, axis=1)
        df_for_clono_pro['chain_cdr3nt'] = df_for_clono_pro.loc[:,['chain', 'cdr3_nt']].apply(':'.join, axis=1)

        cbs = set(df_for_clono_pro['barcode'])
        clonotypes = open(f'{self.outdir}/clonotypes.csv', 'w')
        clonotypes.write('barcode\tcdr3s_aa\tcdr3s_nt\n')
        for cb in cbs:
            temp = df_for_clono_pro[df_for_clono_pro['barcode']==cb]
            temp = temp.sort_values(by='chain', ascending=True)
            aa_chain = ';'.join(list(temp['chain_cdr3aa']))
            nt_chain = ';'.join(list(temp['chain_cdr3nt']))
            clonotypes.write(f'{cb}\t{aa_chain}\t{nt_chain}\n')
        clonotypes.close() 

        df_clonotypes = pd.read_csv(f'{self.outdir}/clonotypes.csv', sep='\t', index_col=None)
        contig_with_clonotype = copy.deepcopy(df_clonotypes)
        df_dict = df_clonotypes[["cdr3s_nt", "cdr3s_aa"]].set_index("cdr3s_nt").to_dict(orient='dict')['cdr3s_aa']
        df_clonotypes = df_clonotypes.groupby('cdr3s_nt', as_index=False).agg({'barcode': 'count'})
        df_clonotypes.rename(columns={'barcode': 'frequency'}, inplace=True)
        sum_f = df_clonotypes['frequency'].sum()

        df_clonotypes['proportion'] = df_clonotypes['frequency'].apply(lambda x: x/sum_f)
        df_clonotypes.sort_values(by='frequency', ascending=False, inplace=True)
        df_clonotypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonotypes.shape[0]+1)]
        df_clonotypes['cdr3s_aa'] = df_clonotypes['cdr3s_nt'].apply(lambda x:df_dict[x])
        df_clonotypes = df_clonotypes.reindex(columns=['clonotype_id', 'frequency', 'proportion', 'cdr3s_aa', 'cdr3s_nt'])
        df_clonotypes.to_csv(f'{self.outdir}/clonotypes.csv', sep=',', index=False) 
        used_for_merge = df_clonotypes[['cdr3s_nt','clonotype_id']]

        df_merge = pd.merge(used_for_merge, contig_with_clonotype, on='cdr3s_nt', how='outer')
        df_merge = df_merge[['barcode', 'clonotype_id']]
        df_all_contig = pd.merge(df_merge, df, on='barcode',how='outer')
        df_all_contig.fillna('',inplace = True)
        df_all_contig = df_all_contig[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
        df_filter_contig = df_all_contig[df_all_contig['barcode'].isin(cell_barcodes)]
        for _df in [df_all_contig, df_filter_contig]:
            _df.loc[~_df.contig_id.isin(filtered_contig), 'clonotype_id'] = ''

        df_all_contig.to_csv(f'{self.outdir}/{self.sample}_all_contig.csv', sep=',', index=False)
        df_filter_contig.to_csv(f'{self.outdir}/{self.sample}_filtered_contig.csv', sep=',', index=False)

    @utils.add_log
    def gen_summary(self, df_for_clono):
        """ Generate metrics in html 
        """
        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        cell_barcodes = set(df_for_clono_pro['barcode'])
        total_cells =len(cell_barcodes)

        read_count = 0
        read_count_all = 0
        umi_dict = defaultdict(set)
        umi_count = defaultdict()
        with pysam.FastxFile(self.fq2) as fq:
            for read in fq:
                read_count_all+=1
                cb = read.name.split('_')[0]
                umi = read.name.split('_')[1]
                umi_dict[cb].add(umi)
                if cb in cell_barcodes:
                    read_count+=1
        for cb in umi_dict:
            umi_count[cb] = len(umi_dict[cb])
        df_umi = pd.DataFrame.from_dict(umi_count, orient='index', columns=['UMI'])
        df_umi['barcode'] = df_umi.index
        df_umi = df_umi.reset_index(drop=True)
        df_umi = df_umi.reindex(columns=['barcode', 'UMI'])
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in cell_barcodes else 'UB')
        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))

        self.add_metric(
            name = 'Estimated Number of Cells',
            value = total_cells,
            help_info = "Number of cells which contain at least one chain (for TCR: TRA or TRB, for BCR: IGH, IGL or IGK)"            
        )

        self.add_metric(
            name = 'Mean Read Pairs per Cell',
            value = int(read_count/total_cells),
            help_info = 'Number of input reads divided by the estimated number of cells'
        )

        used_read = 0
        with pysam.FastxFile(self.assembled_fa) as fa:
            for read in fa:
                bc = read.name.split('_')[0]
                if bc in cell_barcodes:
                    used_read += 1
        self.add_metric(
            name = 'Mean Used Read Pairs per Cell',
            value = int(used_read/total_cells), 
            help_info = "Mean number of reads used in assembly per cell-associated barcode"
        )
        self.add_metric(
            name = 'Fraction of Reads in Cells',
            value = used_read,
            total = read_count_all,
            help_info = 'Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes'
        )
        
        for c in self.chains:
            temp_df = df_for_clono_pro[df_for_clono_pro['chain']==c]

            try:
                median_umi = int(temp_df['umis'].median())
            except ValueError:
                # ValueError: cannot convert float NaN to integer
                median_umi = 0

            self.add_metric(
                name = f'Median {c} UMIs per Cell',
                value = median_umi
            )

    def run(self):
        original_df = self.parse_contig_file()
        df_for_clono, cell_barcodes, filtered_contig = self.cell_calling(original_df)
        self.filter_fasta(cell_barcodes)
        self.parse_clonotypes(original_df, df_for_clono, cell_barcodes, filtered_contig)
        self.gen_summary(df_for_clono)


@utils.add_log
def summarize(args):
    with Summarize(args, display_title="Cells") as runner:
        runner.run()


def get_opts_summarize(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--ref', help='reference name', choices=["hg19", "hg38", "GRCm38", "other"], required=True)
    parser.add_argument('--coef', help='coef for auto filter', default=5)
    parser.add_argument(
        '--diffuseFrac', 
        help="If cell A's two chains CDR3s are identical to another cell B, and A's chain abundance is significantly lower than B's, filter A.",
        action='store_true')
    parser.add_argument(
        "--expected_target_cell_num", 
        help="Expected T or B cell number. If `--target_cell_barcode` is provided, this argument is ignored.", 
        type=int,
        default=3000,
    )
    parser.add_argument(
        '--target_cell_barcode', 
        help="Barcode of target cells. Auto or path of plain text file with one barcode per line",
        default=None)
    parser.add_argument(
        "--target_weight", 
        help="UMIs of the target cells are multiplied by this factor. Only used when `--target_cell_barcode` is provided.", 
        type=float,
        default=6.0,
    )
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq2', help='Barcode R2 reads.', required=True)
        parser.add_argument('--assemble_out', help='Result of  assemble dirctory.', required=True)
        parser.add_argument('--match_dir', help='Match scRNA-seq directory.', required=True)
