from celescope.dynaseq.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_dynaseq(Multi):

    def conversion(self, sample):
        step = 'conversion'
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam'
        cell = f'{self.outdir_dic[sample]["count"]}/{sample}_matrix_10X/barcodes.tsv'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--cell {cell} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def substitution(self, sample):
        step = 'substitution'
        bam = f'{self.outdir_dic[sample]["conversion"]}/{sample}.PosTag.bam'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
        )
        self.process_cmd(cmd, step, sample, m=1, x=1)

    def replacement(self, sample):
        step = 'replacement'
        bam = f'{self.outdir_dic[sample]["conversion"]}/{sample}.PosTag.bam'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--bg {self.col5_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=10, x=1)

    def replace_tsne(self, sample):
        step = 'replace_tsne'
        tsne_file = f'{self.outdir_dic[sample]["analysis"]}/{sample}_tsne_coord.tsv'
        mat_file = f'{self.outdir_dic[sample]["replacement"]}/{sample}.fraction_of_newRNA_matrix.txt'
        rep_file = f'{self.outdir_dic[sample]["replacement"]}/{sample}.fraction_of_newRNA_per_cell.txt'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--tsne {tsne_file} '
            f'--mat {mat_file} '
            f'--rep {rep_file} '
        )
        self.process_cmd(cmd, step, sample, m=1, x=1)


def main():
    multi = Multi_dynaseq(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
