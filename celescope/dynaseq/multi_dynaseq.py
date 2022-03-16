from celescope.dynaseq.__init__ import __ASSAY__
from celescope.tools.multi import Multi
from celescope.tools.__init__ import FILTERED_MATRIX_DIR_SUFFIX, BARCODE_FILE_NAME


class Multi_dynaseq(Multi):

    """
    ## Usage
    ```
        multi_dynaseq\\
        --mapfile ./rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --STAR_param "--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outSAMattributes MD"\\
        --strand /SGRNJ03/Public/Database/genome/gene.strandedness.csv
    ```

    You need to generate strandness-file from gtf file. 
    The format is "geneID,strand", eg:
    ```
    ENSG00000223972,+
    ENSG00000227232,-
    ENSG00000278267,-
    ```
    """


    def conversion(self, sample):
        step = 'conversion'
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam'
        cell = f'{self.outdir_dic[sample]["count"]}/{sample}_{FILTERED_MATRIX_DIR_SUFFIX[0]}/{BARCODE_FILE_NAME[0]}'
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
