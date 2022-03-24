from celescope.tools.multi import Multi
from celescope.convert10X.__init__ import __ASSAY__


class Multi_convert10X(Multi):
    def convert(self, sample):
        step = 'convert'
        cmd_line = self.get_cmd_line(step, sample)
        fq2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--fq2 {fq2} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_convert10X(__ASSAY__)
    multi.run()
    

if __name__ == '__main__':
    main()
