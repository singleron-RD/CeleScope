from celescope.capture_virus.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_capture_virus(Multi):

    def star_virus(self, sample):
        step = 'star_virus'
        cmd_line = self.get_cmd_line(step, sample)
        if self.args.not_consensus:
            fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        else:
            fq = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fq'
            cmd_line += ' --consensus_fq '
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)
    
    def featureCounts_capture_virus(self, sample):
        step = 'featureCounts_capture_virus'
        cmd_line = self.get_cmd_line(step, sample)
        filter_umi_file = f'{self.outdir_dic[sample]["count_capture_virus"]}/{sample}_virus_UMI_count.tsv'
        filter_read_count_json = f'{self.outdir_dic[sample]["count_capture_virus"]}/{sample}_filter_read_count.tsv'
        bam = f'{self.outdir_dic[sample]["star_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--filter_umi_file {filter_umi_file} '
            f'--filter_read_count_json {filter_read_count_json} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)        

    def count_capture_virus(self, sample):
        step = 'count_capture_virus'
        cmd_line = self.get_cmd_line(step, sample)
        virus_bam = f'{self.outdir_dic[sample]["star_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
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
    
    def count_capture_virus_mtx(self, sample):
        step = 'count_capture_virus_mtx'
        cmd_line = self.get_cmd_line(step, sample)       
        if self.args.umi_threshold == "otsu":
            filter_umi_file = f'{self.outdir_dic[sample]["analysis_capture_virus"]}/{sample}_otsu_UMI_count.tsv'
        elif self.args.umi_threshold == "auto":
            filter_umi_file = f'{self.outdir_dic[sample]["analysis_capture_virus"]}/{sample}_auto_UMI_count.tsv'
        bam = f'{self.outdir_dic[sample]["featureCounts_capture_virus"]}/{sample}_filter_name_sorted.bam'
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--filter_umi_file {filter_umi_file} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample)    


def main():
    multi = Multi_capture_virus(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
