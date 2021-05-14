from celescope.vdj.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_vdj(Multi):

    def mapping_vdj(self, sample):
        step = 'mapping_vdj'
        cmd_line = self.get_cmd_line(step, sample)
        if self.args.not_consensus:
            fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        else:
            fq = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fq'
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)


    def count_vdj(self, sample):
        # count_vdj
        step = 'count_vdj'
        cmd_line = self.get_cmd_line(step, sample)
        UMI_count_filter1_file = (
            f'{self.outdir_dic[sample]["mapping_vdj"]}/{sample}_UMI_count_filtered1.tsv'
        )
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--UMI_count_filter1_file {UMI_count_filter1_file} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)



def main():
    multi = Multi_vdj(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()

