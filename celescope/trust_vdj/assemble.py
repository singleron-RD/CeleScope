import os
from celescope.tools import utils
from celescope.tools.Step import Step, s_common
import pysam
import pandas as pd
from collections import defaultdict
import glob
import re
from Bio.Seq import Seq


TRUST = '/SGRNJ03/randd/zhouxin/software/TRUST4/'

@utils.add_log
def count_fq(fq):
    dic = defaultdict(list)
    with pysam.FastxFile(fq) as fq:
        for entry in fq:
            attr = entry.sequence
            cb = attr[:24]
            umi = attr[24:]
            name = entry.name
            dic['barcode'].append(cb)
            dic['UMI'].append(umi)
            dic['seq_name'].append(name)

    count_df = pd.DataFrame(dic, columns=list(dic.keys()))

    return count_df


class Assemble(Step):
    """
    Features

    - Get fq file
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.sample = args.sample
        self.species = args.species
        self.speed_up = args.speed_up
        self.match_dir = args.match_dir
        self.cells = args.cells

    @utils.add_log
    def get_barcodes(self):        
        tsne = glob.glob(f'{self.match_dir}/06.analysis/*_tsne_coord.tsv')
        tsne = tsne[0]
        tsne_coord = pd.read_csv(tsne, sep='\t', index_col=0)
        barcodes = tsne_coord.index.tolist()

        # write barcodes
        res = [] 
        for barcode in barcodes:
            barcode = Seq(barcode)
            barcode_reversed = barcode.reverse_complement()
            bc = str(barcode_reversed)
            res.append(bc)

        df = pd.DataFrame(res, columns=['barcode'])

        return df

    @utils.add_log
    def cut_off(self):
        barcodes = self.get_barcodes()
        df = count_fq(self.fq1)
        df_umi = df.groupby(['barcode', 'UMI'], as_index=False).agg({'seq_name': 'count'})
        df_umi = df_umi.groupby(['barcode'], as_index=False).agg({'UMI': 'count'})

        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi = df_umi.reset_index()

        UMI_num = int(self.cells)
        rank = UMI_num / 100
        rank_UMI = df_umi.loc[rank, 'UMI']
        UMI_min = int(rank_UMI / 10)

        df_umi_filtered = df_umi[df_umi.UMI >= UMI_min]

        df_tmp = pd.merge(df_umi_filtered, barcodes, on='barcode', how='inner')

        matched_barcodes = df_tmp.barcode.tolist()
        with open(f'{self.outdir}/{self.sample}_matched_barcodes.txt', 'w') as fh:
            for barcode in matched_barcodes:
                fh.write(str(barcode)+ '\n')
        string = f'Get {len(matched_barcodes)} matched barcodes'

        Assemble.cut_off.logger.info(string)

        df_all = pd.merge(df_tmp, df, on='barcode', how='outer')
        seq_list = df_all['seq_name'].tolist()

        with open(f'{self.outdir}/seqlist.txt', 'w') as fh:
            for name in seq_list:
                fh.write(str(name)+'\n')

    
    @utils.add_log
    def getFqfile(self):

        cmd1 = (
            f'seqtk subseq {self.fq1} {self.outdir}/seqlist.txt > {self.outdir}/{self.sample}_matched_R1.fq'
        )
        os.system(cmd1)

        cmd2 = (
            f'seqtk subseq {self.fq2} {self.outdir}/seqlist.txt > {self.outdir}/{self.sample}_matched_R2.fq'
        )
        os.system(cmd2)


    @utils.add_log
    def run(self):

        self.cut_off()
        self.getFqfile()

        species = self.species

        index_file = f'{TRUST}/index/{species}/{species}_ref.fa'
        ref = f'{TRUST}/index/{species}/{species}_IMGT+C.fa'

        string1 = ''
        if self.speed_up:
            string1 = '--repseq '
        cmd = (
            f'{TRUST}/run-trust4 -t {self.thread} '
            f'-u {self.outdir}/{self.sample}_matched_R2.fq '
            f'--barcode {self.outdir}/{self.sample}_matched_R1.fq '
            f'--barcodeRange 0 23 + '
            f'-f {index_file} '
            f'--ref {ref} '
            f'{string1}'
            f'-o {self.sample} --od {self.outdir}/TRUST4' 
        )

        Assemble.run.logger.info(cmd)

        if not os.path.exists(f'{self.outdir}/TRUST4/{self.sample}_barcode_report.tsv'):
            os.system(cmd)

            #fq = f'{self.outdir}/TRUST4/{self.sample}_toassemble.fq'

        # report
        os.system(f'rm {self.outdir}/seqlist.txt')


@utils.add_log
def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='R1 reads from barcode step', required=True)
        parser.add_argument('--fq2', help='R2 reads from barcode step', required=True)
        parser.add_argument('--match_dir', help='rna analysis dir', required=True)

    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)
    parser.add_argument('--cells', help='expected cell number', default=3000)
    parser.add_argument('--speed_up', help='speed assemble for TCR/BCR seq data', action='store_true')       








