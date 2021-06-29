from logging import raiseExceptions
import os
from celescope.tools import utils
from celescope.tools.step import Step, s_common
import pysam
import pandas as pd
from collections import defaultdict
import glob
from Bio.Seq import Seq


@utils.add_log
def count_fq(fq):
    dic = defaultdict(list)
    with pysam.FastxFile(fq) as fq:
        for entry in fq:
            attr = entry.name.split('_')
            cb = attr[0]
            umi = attr[1]
            name = entry.name
            dic['barcode'].append(cb)
            dic['UMI'].append(umi)
            dic['seq_name'].append(name)

        count_df = pd.DataFrame(dic, columns=list(dic.keys()))

    return count_df


class Matching(Step):
    """
    Features

    - Cut off V(D)J data by UMI count. Default value is 1/10 of the 30th barcode's UMIs ranked by UMI count.
    - Match V(D)J barcodes after cut off with RNA cell barcodes.

    Output

    - `02.matching/count.txt`. Record the UMI count of each barcode in raw V(D)J data.
    - `02.matching/{sample}_matched_barcodes.txt`. Contain the matched barcode.
    - `02.matching/{sample}_matched_R1.fq`, `02.match/{sample}_matched_R2.fq. Barcode and UMI are contained in the R1 reads.

    """
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.sample = args.sample
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
        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)

        UMI_num = int(self.cells)
        rank = UMI_num / 100
        rank_UMI = df_umi.loc[rank, 'UMI']
        UMI_min = int(rank_UMI / 10)

        df_umi_filtered = df_umi[df_umi.UMI >= UMI_min]

        df_tmp = pd.merge(df_umi_filtered, barcodes, on='barcode', how='inner')

        #output rna barcodes
        matched_barcodes = df_tmp.barcode.tolist()
        with open(f'{self.outdir}/{self.sample}_matched_barcodes.txt', 'w') as fh:
            for barcode in matched_barcodes:
                cb = Seq(barcode)
                cb = cb.reverse_complement()
                cb = str(cb)
                fh.write(cb+ '\n')
        string = f'Get {len(matched_barcodes)} matched barcodes'

        Matching.cut_off.logger.info(string)

        return matched_barcodes

    
    @utils.add_log
    def getFqfile(self):

        bcs = self.cut_off()
        matched_fq = open(f'{self.outdir}/{self.sample}.bam', 'w')
        if len(bcs) > 0:
            with pysam.AlignmentFile(self.bamfile) as bam:
                for read in bam:
                    attrs = read.reference_name.split('_')
                    cb = attrs[0]
                    if cb in bcs:
                        matched_fq.write(str(read)+'\n')
        else:
            raise Exception('No matched barcodes!')


    @utils.add_log
    def run(self):
        self.getFqfile()


@utils.add_log
def matching(args):
    step_name = 'matching'
    match_obj = Matching(args, step_name)
    match_obj.run()


def get_opts_matching(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--match_dir', help='rna analysis dir', required=True)
        parser.add_argument('--fq1', help='R1 reads from convert step', required=True)
        parser.add_argument('--fq2', help='R2 reads from convert step', required=True)
    parser.add_argument('--cells', help='expected cell number', default=3000)


