from celescope.tools.multi import Multi
from celescope.vdj10X.__init__ import __ASSAY__


class Multi_vdj10X(Multi):

    def convert(self, sample):
        step = 'convert'
        cmd_line = self.get_cmd_line(step, sample)
        fq2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--fq2 {fq2} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fqs_dir = f'{self.outdir_dic[sample]["convert"]}'
        barcode_dic = f'{fqs_dir}/barcode_correspond.txt'
        match_dir = f'{self.col4_dict[sample]}'
        cmd = (
            f'{cmd_line} '
            f'--fqs_dir {fqs_dir} '
            f'--match_dir {match_dir} '
            f'--barcode_dic {barcode_dic} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.mem, x=self.args.thread)
        
    def match(self, sample):
        step = 'match'
        cmd_line = self.get_cmd_line(step, sample)
        match_dir = f'{self.col4_dict[sample]}'
        contig_df = f'{self.outdir_dic[sample]["assemble"]}/contigs.csv'
        cmd = (
            f'{cmd_line} '
            f'--match_dir {match_dir} '
            f'--contig_df {contig_df} ' 
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

def main():
    multi = Multi_vdj10X(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()

