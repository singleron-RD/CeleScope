from celescope.capture_rna.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_capture_rna(Multi):

    def count_capture_rna(self, sample):
        step = 'count_capture_rna'
        cmd_line = self.get_cmd_line(step, sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_name_sorted.bam'
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=10, x=1)

    def analysis(self, sample):
        step = 'analysis'
        cmd_line = self.get_cmd_line(step, sample)
        out_dir = f'{self.outdir_dic[sample][step]}'
        matrix_file = f'{self.outdir_dic[sample]["count_capture_rna"]}/{sample}_matrix.tsv.gz'
        cmd = (
            f'{cmd_line} '
            f'--matrix_file {matrix_file} '
        )
        self.process_snakemake_cmd(cmd, step, out_dir,sample,x=1)
        self.process_cmd(cmd, step, sample, m=10, x=1)


def main():
    multi = Multi_capture_rna(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
