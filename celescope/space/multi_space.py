from celescope.space.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_space(Multi):
    def starsolo(self, sample):
        step = "starsolo"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]}  --spatial {self.col4_dict[sample]} '
        )
        self.process_cmd(
            cmd,
            step,
            sample,
            m=int(self.args.limitBAMsortRAM / 1e9),
            x=self.args.thread,
        )

    def analysis(self, sample):
        step = "analysis"
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f"{cmd_line} "
            f"--raw {self.outdir_dic[sample]['outs']}/raw "
            f"--filtered {self.outdir_dic[sample]['outs']}/filtered "
            f"--spatial {self.col4_dict[sample]} "
        )
        self.process_cmd(cmd, step, sample)


def main():
    multi = Multi_space(__ASSAY__, min_col=4, pair=(1, 2))
    multi.run()


if __name__ == "__main__":
    main()
