from celescope.vdj10X.__init__ import  __ASSAY__
from celescope.tools.Multi import Multi
from celescope.tools.__init__ import __PATTERN_DICT__


class Multi_vdj10X(Multi):

    def convert(self, sample):
        step = 'convert'
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


    def vdj_10X(self, sample):
        step = 'vdj_10X'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.mem, x=self.args.thread)



def main():
    multi = Multi_vdj10X(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()

