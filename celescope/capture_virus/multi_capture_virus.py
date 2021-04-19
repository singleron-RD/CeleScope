from celescope.capture_virus.__init__ import __ASSAY__
from celescope.tools.Multi import Multi


class Multi_capture_virus(Multi):

    def STAR_virus(self, sample):
        step = 'STAR_virus'
        cmd_line = self.get_cmd_line(step, sample)
        input_read = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{cmd_line} '
            f'--input_read {input_read} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def count_capture_virus(self, sample):
        step = 'count_capture_virus'
        cmd_line = self.get_cmd_line(step, sample)
        virus_bam = f'{self.outdir_dic[sample]["STAR_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{cmd_line} '
            f'--virus_bam {virus_bam} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def analysis_capture_virus(self, sample):        
        step = 'analysis_capture_virus'
        cmd_line = self.get_cmd_line(step, sample)
        virus_file = f'{self.outdir_dic[sample]["count_capture_virus"]}/{sample}_virus_UMI_count.tsv'
        cmd = (
            f'{cmd_line} '
            f'--virus_file {virus_file} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_capture_virus(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()



