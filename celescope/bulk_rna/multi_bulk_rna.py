from celescope.bulk_rna.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_bulk_rna(Multi):
    def starsolo(self, sample):
        step = "starsolo"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)


def main():
    multi = Multi_bulk_rna(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
