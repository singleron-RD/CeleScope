from celescope.tools.multi import Multi
from celescope.vdj_full_len.__init__ import __ASSAY__

class Multi_vdj_full_len(Multi):

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
    
    def check(self, sample):
        step = 'check'
        cmd_line = self.get_cmd_line(step, sample)
        rep = f'{self.outdir_dic[sample]["assemble"]}/stat.txt'
        rep_bc = f'{self.outdir_dic[sample]["barcode"]}/stat.txt'
        html = f'{sample}/{sample}_report.html'
        cmd = (
            f'{cmd_line} '
            f'--rep {rep} '
            f'--rep_bc {rep_bc} '
            f'--html {html}'
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def mapping(self,sample):
        step = 'mapping'
        cmd_line = self.get_cmd_line(step,sample)
        match_dir = f'{self.col4_dict[sample]}'
        cmd=(
            f'{cmd_line} '
            f'--match_dir {match_dir} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

def main():
    multi = Multi_vdj_full_len(__ASSAY__)
    multi.run()
    

if __name__ == '__main__':
    main()

