import os,re
import sys
import pandas as pd
import anndata
import pysam
from multiprocessing import Pool
import subprocess

from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.dynaseq.__init__ import DYNA_MATRIX_DIR_SUFFIX
from celescope.rna.mkref import Mkref_rna
from celescope.tools.matrix import CountMatrix, Features, ROW, COLUMN
from celescope.tools import reference
from celescope.tools.plotly_plot import Tsne_plot, Violin_plot

toolsdir = os.path.dirname(__file__)


class Replacement(Step):
    """
    Features
    - Quantify unlabeled and labeled RNA.
    - Boxplots for TOR rates distribution.
    - TSNE plot for TOR rate 

    Output
    - `{sample}.labeled.h5ad` h5ad file contains ['total', 'labeled', 'unlabeled'] layers and TOR rate of each cell/gene.
    - `{sample}_labeled_feature_bc_matrix` The labeled expression matrix of cell barcodes & all features in Matrix Market Exchange Formats. (will be deprecated in future versions)
    - `{sample}_unlabeled_feature_bc_matrix` The unlabeled expression matrix of cell barcodes & all features in Matrix Market Exchange Formats. (will be deprecated in future versions)
    - `{sample}_labeled_detail.txt`  tab-delimited  file:
        - Barcode: Cell barcode sequence
        - UMI: UMI sequence
        - geneID: gene ID
        - TC: TC site number in a read (backgroup snp removed)
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # input files
        self.outdir = args.outdir
        self.sample = args.sample
        self.thread = int(args.thread)
        self.inbam = args.bam
        self.bcfile = args.cell
        self.snp_file = args.bg
        self.tsne = args.tsne
        self.cellsplit = args.cellsplit

        # set
        self.cell_dict, self.cell_num = utils.barcode_list_stamp(self.bcfile,cut=self.cellsplit)
        gtf_file = Mkref_rna.parse_genomeDir(args.genomeDir)['gtf']
        gp = reference.GtfParser(gtf_file)
        gp.get_id_name()
        self.features = gp.get_features()
        self.used_gene_id, self.used_gene_name = [], []
        self.used_features = None
        self.totaldf = pd.DataFrame()
        self.newdf, self.olddf = pd.DataFrame(), pd.DataFrame()
        self.adata = anndata.AnnData()
        self.cell_list = []
        self.bg = None

        # output files
        ## tmp outputdir
        self.tmp_dir = f'{args.outdir}/tmp'
        utils.check_mkdir(self.tmp_dir)
        ## final outputs
        self.h5ad = f'{self.out_prefix}.labeled.h5ad'
        self.detail_txt = f'{self.out_prefix}_labeled_detail.txt'
        self.dir_labeled = f'{self.out_prefix}_{DYNA_MATRIX_DIR_SUFFIX[0]}'
        self.dir_unlabeled = f'{self.out_prefix}_{DYNA_MATRIX_DIR_SUFFIX[1]}'


    @utils.add_log
    def run(self):
        # get backgroud snp        
        self.bg = self.background_snp()
        # replacement
        self.run_quant()
        # stat and plot
        self.tor_plot()
        self.add_help()
        # output dedup and clean 
        self.clean_tmp()


    @utils.add_log
    def run_quant(self):
        ## set Parallelism para
        cell_arr = []
        fetch_arr = [self.inbam] * len(self.cell_dict)
        snp_list = [self.bg] * len(self.cell_dict)
        for x in self.cell_dict:
            cell_arr.append(self.cell_dict[x])
            self.cell_list.extend(self.cell_dict[x])
    
        mincpu = min(self.cell_num, self.thread)
        with Pool(mincpu) as pool:
            results = pool.starmap(Replacement.modify_bam, 
                      zip(fetch_arr,snp_list,cell_arr))
        # merge matrix
        for i in results:
            self.totaldf = pd.concat([self.totaldf,i])
        self.newdf = self.totaldf[self.totaldf['TC']>0]
        self.olddf = self.totaldf[self.totaldf['TC']==0]
        self.totaldf.to_csv(self.detail_txt, sep="\t", index=False)

        used_gene = self.totaldf['geneID'].unique()        
        for index,item in enumerate(self.features.gene_id):
            if item in used_gene:
                self.used_gene_id.append(item)
                self.used_gene_name.append(self.features.gene_name[index])
        self.used_features = Features(self.used_gene_id,self.used_gene_name)

        # output
        tmp_newdf = self.newdf.groupby([COLUMN, ROW]).agg({'UMI': 'count'})
        tmp_olddf = self.olddf.groupby([COLUMN, ROW]).agg({'UMI': 'count'})
        self.write_sparse_matrix(tmp_newdf, self.dir_labeled)
        self.write_sparse_matrix(tmp_olddf, self.dir_unlabeled)
        self.write_h5ad(self.totaldf, self.newdf, self.olddf)


    @staticmethod
    def modify_bam(bam, bg, cells):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, 'rb')
        pysam.set_verbosity(save)
        readdict= {}

        for read in bamfile.fetch(until_eof=True):
            try:
                cb = read.get_tag('CB')
                if cb not in cells:
                    continue
                chro = read.reference_name
                ub = read.get_tag('UB')
                gene = read.get_tag('GX')
                tctag = 0
                true_tc = []

                if read.get_tag('ST') == '+':
                    stag = read.get_tag('TL')
                else:
                    stag = read.get_tag('AL')
                if len(stag) == 1 and stag[0] == 0:
                    tctag = 0
                    true_tc = stag
                else:
                    for si in range(0, len(stag)):
                        pos = chro + '_' + str(stag[si])
                        if pos not in bg:
                            true_tc.append(int(stag[si]))
                    tctag = len(true_tc)
                ## dedup: select the most TC read per UMI_gene
                readid = f'{cb}_{ub}_{gene}'
                if readid not in readdict:
                    readdict[readid] = [tctag,read,true_tc]
                else:
                    if tctag > readdict[readid][0]:
                        readdict[readid] = [tctag,read,true_tc]

            except (ValueError, KeyError):
                continue
        bamfile.close()

        ## count df
        tc_df = pd.DataFrame.from_dict(readdict,orient='index', columns=["TC","read","loc"])
        tc_df.reset_index(inplace=True)
        ub_df = tc_df['index'].str.split('_',expand=True)
        ub_df.columns = ['Barcode','UMI','geneID']
        #ub_df['Barcode'] = [cell] * ub_df.shape[0]
        outframe = pd.concat([ub_df,tc_df['TC']], axis=1)

        return outframe

    @staticmethod
    def createTag(d):
        return ''.join([''.join(key) + str(d[key]) + ';' for key in d.keys()])[:-1]

    @staticmethod
    def modifySCTag(sc,cnt,sstag):
        pattern = re.compile(rf'{sstag}\d+')
        result = pattern.sub(f'{sstag}{cnt}', sc)
        return result

    @utils.add_log
    def background_snp(self):
        outdict = {}
        bgs = []
        for bgargv in self.snp_file:
            if ',' in bgargv:
                bgs += bgargv.strip().split(',')
            else:
                bgs.append(bgargv)
        
        for bgfile in bgs:
            if bgfile.endswith('.csv'):
                df = pd.read_csv(bgfile, dtype={"chrom":str})
                if 'pos' in df.columns:
                    df['chrpos'] = df['chrom']+'_'+df['pos'].astype(str)
                else: #compatible with previous version
                    df['chrpos'] = df['chrom']+'_'+df['pos2'].astype(str)
                df1 = df[['chrpos','convs']]
                df1.set_index('chrpos',inplace=True)
                for key1 in df1.index.to_list():
                    outdict[key1] = 1
                #outdict.update(df1.to_dict(orient='index'))
            elif bgfile.endswith('.vcf'):
                bcf_in = pysam.VariantFile(bgfile)
                for rec in bcf_in.fetch():
                    try:
                        chrom, pos = rec.chrom, rec.pos
                        chr_pos = chrom+'_'+str(pos-1)
                        outdict[chr_pos] = 1
                    except (ValueError, KeyError):
                        continue
                bcf_in.close()
            else:
                try:
                    sys.exit(1)
                except SystemExit:
                    print('Background snp file format cannot be recognized! Only csv or vcf format.')
                finally:
                    print('Background snp file format cannot be recognized! Only csv or vcf format.')
        return outdict


    @utils.add_log
    def write_sparse_matrix(self, df, matrix_dir):
        count_matrix = CountMatrix.from_dataframe(df, self.features, barcodes=self.cell_list)
        count_matrix.to_matrix_dir(matrix_dir)

    @utils.add_log
    def write_h5ad(self, df, df_new, df_old):
        layers = {}
        matrix = CountMatrix.dataframe_to_matrix(df, self.used_features, barcodes=self.cell_list)
        layers['total'] = matrix
        layers['labeled'] = CountMatrix.dataframe_to_matrix(df_new, self.used_features, barcodes=self.cell_list)
        layers['unlabeled'] = CountMatrix.dataframe_to_matrix(df_old, self.used_features, barcodes=self.cell_list)
        
        adata = anndata.AnnData(
            X=matrix,
            obs=pd.DataFrame(index=pd.Series(self.cell_list, name='cell')),
            var=pd.DataFrame(index=pd.Series(self.used_gene_id, name='gene_id'), data={'gene_name': pd.Categorical(self.used_gene_name)}),
            layers=layers
        )
        adata = self.tor_stat(adata)
        adata.write(self.h5ad)
        self.adata = adata

    @utils.add_log
    def tor_stat(self, adata):
        cell_ntr = adata.layers['labeled'].sum(axis=1) / adata.layers['total'].sum(axis=1)
        gene_ntr = adata.layers['labeled'].sum(axis=0) / adata.layers['total'].sum(axis=0)
        adata.obs['TOR'] = cell_ntr
        adata.var['TOR'] = gene_ntr.T
        return adata

    @utils.add_log
    def tor_plot(self):
        tsne = pd.read_csv(self.tsne,sep="\t",index_col=0)
        tsne = tsne.loc[self.adata.obs.index]
        self.adata.obs['tSNE_1'] = tsne['tSNE_1']
        self.adata.obs['tSNE_2'] = tsne['tSNE_2']

        tsne_tor = Tsne_plot(self.adata.obs.sort_values(by="TOR"), 'TOR', discrete=False)
        tsne_tor.set_color_scale("PuRd")        
        self.add_data(tsne_tor=tsne_tor.get_plotly_div())

        vln_gene = Violin_plot(self.adata.var['TOR'], 'gene', color='#1f77b4').get_plotly_div()
        self.add_data(violin_gene=vln_gene)
        vln_cell = Violin_plot(self.adata.obs['TOR'], 'cell', color='#ff7f0e').get_plotly_div()
        self.add_data(violin_cell=vln_cell)

    @utils.add_log
    def add_help(self):
        self.add_help_content(
            name='TOR:',
            content='(RNA turn-over rate) Fraction of labeled transcripts per gene or cell.'
        )


    def run_cmd(self, cmd):
        subprocess.call(' '.join(cmd), shell=True)

    @utils.add_log
    def clean_tmp(self):
        cmd = (f'rm -rf {self.tmp_dir}')
        self.debug_subprocess_call(cmd)


@utils.add_log
def replacement(args):
    if args.control:
        return
    with Replacement(args, display_title='Labeled') as runner:
        runner.run()


def get_opts_replacement(parser, sub_program):
    parser.add_argument('--genomeDir',help=HELP_DICT['genomeDir'])
    parser.add_argument("--control", action='store_true',
                        help="For control samples to generate backgroup snp files and skip replacement")
    parser.add_argument("--cellsplit", default=300, type=int, help='split N cells into a list')
    if sub_program:
        parser.add_argument("--bam", 
            help='convsrion tagged bam from conversion step', required=True)
        parser.add_argument("--cell", help='barcode cell list', required=True)
        parser.add_argument('--bg', nargs='+', required=False,
                            help='background snp file, csv or vcf format')
        parser.add_argument('--tsne', help='tsne file from analysis step', required=True)
        parser = s_common(parser)
    return parser
