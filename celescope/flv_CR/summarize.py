import json
import pandas as pd
import pysam
import os 

from celescope.tools import utils
from celescope.tools.step import s_common, Step


class Summarize(Step):
    """
    ## Features

    - Convert 10X barcode of assemble result back to SGR barcode.

    - Generate Productive contigs sequences and annotation files.

    - Generate VDJ-annotation metrics in html.

    ## Output

    - `filtered_contig_annotations.csv` High-level annotations of each high-confidence contigs from cell-associated barcodes.

    - `filtered_contig.fasta` High-confidence contig sequences annotated in the filtered_contig_annotations.csv.

    - `productive_contig_annotations.csv` Annotations of each productive contigs from cell-associated barcodes. This is a subset of filtered_contig_annotations.csv.

    - `productive_contig.fasta` Productive contig sequences annotated in the productive_contig_annotations.csv.

    - `clonotypes.csv` High-level descriptions of each clonotype.

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.version = os.path.dirname(args.soft_path).split('/')[-1].split('-')[-1]

        with open(args.barcode_convert_json, 'r') as f:
            self.tenX_sgr = json.load(f)

        # assemble result directory in 10X barcode.
        self.assemble_out = args.assemble_out

        # out
        self.filtered_contig_annotations = f'{self.outdir}/filtered_contig_annotations.csv'
        self.filter_contig_fasta = f'{self.outdir}/filtered_contig.fasta'
    
    @utils.add_log
    def convert_barcode_to_SGR(self):
        """
        Convert 10X barcode to SGR barcode format.
        """

        annotation_file = pd.read_csv(f'{self.assemble_out}/filtered_contig_annotations.csv', sep=',', index_col=None)
        annotation_file['barcode'] = annotation_file['barcode'].apply(lambda x: utils.reverse_complement(self.tenX_sgr[x.split('-')[0]]))
        annotation_file['contig_id'] = annotation_file['contig_id'].apply(lambda x: utils.reverse_complement(
            self.tenX_sgr[x.split('-')[0]])+'_'+x.split('_')[1]+'_'+x.split('_')[2])
        if int(self.version.split('.')[0]) < 4:
            annotation_file['productive'].replace({'True': True, 'None': False}, inplace=True)

        annotation_file.to_csv(self.filtered_contig_annotations, sep=',', index=False)
    
        tenX_fasta_file = pysam.FastxFile(f'{self.assemble_out}/filtered_contig.fasta')
        SGR_fasta_file = open(self.filter_contig_fasta, 'w')
        for entry in tenX_fasta_file:
            name = entry.name
            seq = entry.sequence
            attrs = name.split('_')
            new_name = utils.reverse_complement(self.tenX_sgr[attrs[0].split('-')[0]]) + '_' + attrs[1] + '_' + attrs[2]
            SGR_fasta_file.write(f'>{new_name}\n{seq}\n')
        SGR_fasta_file.close()

        os.system(f'cp {self.assemble_out}/clonotypes.csv {self.outdir}')

        Summarize.gen_productive_contig(annotation_file, self.filter_contig_fasta, self.outdir)

    @staticmethod
    @utils.add_log
    def gen_productive_contig(annotation_file, fasta_file, outdir, prefix=''):
        """
        Generate productive_contig_annotation.csv and productive_contig.fasta file.
        """
        productive_contig = annotation_file[annotation_file['productive'] == True]
        productive_contig_id = set(productive_contig['contig_id'])
        productive_contig.to_csv(f'{outdir}/{prefix}productive_contig_annotations.csv', sep=',', index=False)

        productive_fasta = open(f'{outdir}/{prefix}productive_contig.fasta', 'w')
        with pysam.FastxFile(fasta_file, 'r') as f:
            for read in f:
                if read.name in productive_contig_id:
                    productive_fasta.write(">" + read.name + "\n" + read.sequence + "\n")
        productive_fasta.close()

    def run(self):
        self.convert_barcode_to_SGR()


def summarize(args):
    with Summarize(args) as runner:
        runner.run()


def get_opts_summarize(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--soft_path', help='soft path for cellranger')
    if sub_program:
        s_common(parser)
        parser.add_argument('--barcode_convert_json', help='json file', required=True)
        parser.add_argument('--assemble_out', help='directory of cellranger assemble result', required=True)
    return parser