from celescope.ffpe.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_ffpe(Multi):
    def mirna(self, sample):
        step = "mirna"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
        self.process_cmd(cmd, step, sample, m=10, x=self.args.thread)


def main():
    multi = Multi_ffpe(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
