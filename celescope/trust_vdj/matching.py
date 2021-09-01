from collections import defaultdict

import pandas as pd
import pysam
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common


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
        self.fq2 = args.fq2
        self.sample = args.sample
        self.match_dir = args.match_dir

        self.match_barcodes, cell_num = utils.read_barcode_file(self.match_dir)

    @utils.add_log
    def get_match_fastq(self):
        out_fq1 = open(f'{self.outdir}/{self.sample}_matched_R1.fq', 'w')
        out_fq2 = open(f'{self.outdir}/{self.sample}_matched_R2.fq', 'w')
        
        with pysam.FastxFile(self.fq2) as fa:
            for read in fa:
                attr = read.name.split('_')
                cb = attr[0]
                umi = attr[1]
                qual = 'F' * len(cb + umi)
                reversed_cb = str(Seq(cb).reverse_complement())
                if reversed_cb in self.match_barcodes:
                    seq1 = f'@{read.name}\n{cb}{umi}\n+\n{qual}\n'
                    out_fq1.write(seq1)
                    out_fq2.write(str(read)+'\n')

            out_fq1.close()
            out_fq2.close()


    @utils.add_log
    def run(self):
        self.get_match_fastq()


@utils.add_log
def matching(args):
    step_name = 'matching'
    match_obj = Matching(args, step_name)
    match_obj.run()


def get_opts_matching(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--match_dir', help='rna analysis dir', required=True)
        parser.add_argument('--fq2', help='R2 reads from convert step', required=True)


