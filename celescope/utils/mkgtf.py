import sys
import csv
import re
import collections

from celescope.tools import utils


ATTRIBUTES = {
    'gene_biotype':{
        'protein_coding', 
        'lncRNA',
        'antisense',
        'IG_LV_gene',
        'IG_V_gene',
        'IG_V_pseudogene',
        'IG_D_gene',
        'IG_J_gene',
        'IG_J_pseudogene',
        'IG_C_gene',
        'IG_C_pseudogene',
        'TR_V_gene',
        'TR_V_pseudogene',
        'TR_D_gene',
        'TR_J_gene',
        'TR_J_pseudogene',
        'TR_C_gene',
    }
}


class Mkgtf:
    """
    ##Features

    - Filter GTF files.
    
    GTF files downloaded from sites like ENSEMBL and UCSC often contain transcripts and genes which need to be filtered.
    Celescope provides mkgtf, a simple utility to filter genes. The command syntax requires input and output fitlered GTF file names.

    The following gene biotypes will be kept.
    --attribute=gene_biotype:protein_coding 
    --attribute=gene_biotype:lncRNA 
    --attribute=gene_biotype:antisense 
    --attribute=gene_biotype:IG_LV_gene 
    --attribute=gene_biotype:IG_V_gene 
    --attribute=gene_biotype:IG_V_pseudogene 
    --attribute=gene_biotype:IG_D_gene 
    --attribute=gene_biotype:IG_J_gene 
    --attribute=gene_biotype:IG_J_pseudogene 
    --attribute=gene_biotype:IG_C_gene 
    --attribute=gene_biotype:IG_C_pseudogene 
    --attribute=gene_biotype:TR_V_gene 
    --attribute=gene_biotype:TR_V_pseudogene 
    --attribute=gene_biotype:TR_D_gene 
    --attribute=gene_biotype:TR_J_gene 
    --attribute=gene_biotype:TR_J_pseudogene 
    --attribute=gene_biotype:TR_C_gene
    
    """

    def __init__(self, args):

        # in
        self.gtf = args.gtf
        
        # out
        self.filtered_gtf = args.filtered_gtf

    @staticmethod
    def get_properties_dict(properties_str):
        
        if isinstance(properties_str, dict):
            return properties_str

        properties = collections.OrderedDict()
        pattern = re.compile(r'(\S+?)\s*"(.*?)"')
        for m in re.finditer(pattern, properties_str):
            key = m.group(1)
            value = m.group(2)
            properties[key] = value

        return properties

    @staticmethod
    def gtf_reader_iter(raw_gtf):
        with open(raw_gtf, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(reader):
                if len(row) == 0:
                    continue
                if row[0].startswith('#'):
                    yield row, True, None
                    continue

                if len(row) != 9:
                    sys.exit(f"Invalid number of columns in GTF line {i+1}: {row}\n")

                if row[6] not in ['+', '-']:
                    sys.exit(f"Invalid strand in GTF line {i+1}: {row}\n")

                properties = Mkgtf.get_properties_dict(row[8])

                yield row, False, properties


    @utils.add_log
    def mk_gtf(self):
        print("Writing new genes GTF file (may take 10 minutes for a 1GB input GTF file)...")
        with open(self.filtered_gtf, 'w') as f:
            writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
            for row, is_comment, properties in Mkgtf.gtf_reader_iter(self.gtf):
                if is_comment:
                    writer.writerow(row)
                    continue

                remove = False
                for key, value in properties.items():
                    if key in ATTRIBUTES and value not in ATTRIBUTES[key]:
                        remove = True

                if not remove:
                    writer.writerow(row)

            print("...done\n")

    def run(self):
        self.mk_gtf()


def mkgtf(args):
    runner = Mkgtf(args)
    runner.run()


def get_opts_mkgtf(parser, sub_program=True):
    parser.add_argument('gtf', help='raw gtf file')
    parser.add_argument('filtered_gtf', help='filtered gtf file')

    return parser