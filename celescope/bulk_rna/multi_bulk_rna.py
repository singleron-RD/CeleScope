from celescope.bulk_rna.__init__ import __ASSAY__
from celescope.tools.multi import Multi
from celescope.tools.__init__ import (
    OUTS_DIR,
)


class Multi_bulk_rna(Multi):
    def starsolo(self, sample):
        step = "starsolo"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
        self.process_cmd(
            cmd,
            step,
            sample,
            m=int(self.args.limitBAMsortRAM / 1e9),
            x=self.args.thread,
        )

    def split(self, sample):
        step = "split"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr["fq1_str"]} '
            f'--bam {self.outdir_dic[sample][OUTS_DIR]}/{sample}_Aligned.sortedByCoord.out.bam '
            f'--barcodes {self.outdir_dic[sample][OUTS_DIR]}/filtered/barcodes.tsv.gz '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_bulk_rna(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
