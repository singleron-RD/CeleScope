import pysam
from celescope.tools import utils
from celescope.tools.step import Step, s_common


class Convert(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        
        self.fq2 = args.fq2
        
    @utils.add_log  
    def get_fq1(self):
        fq2 = pysam.FastxFile(self.fq2, 'r')
        fq1 = open(f'{self.outdir}/{self.sample}_clean_1.fq', 'w')
        qual = 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
        for entry in fq2:
            name = entry.name
            attrs = name.split('_')
            barcode = attrs[0]
            umi = attrs[1]
            fq1.write(f'@{name}\n{barcode}{umi}\n+\n{qual}\n')
            fq1.flush()
        fq1.close()
        
    def run(self):
        self.get_fq1()
        
def convert(args):
    step_name = 'convert'
    convert_obj = Convert(args, step_name)
    convert_obj.run()
    
def get_opts_convert(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq2', help='clean fq2 file', required=True)
