from celescope.sweetseq.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_sweetseq(Multi):

    def mapping(self, sample):
        step = 'mapping'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '

        )
        self.process_cmd(cmd, step, sample, m=2, x=1)

    def count(self, sample):
        step = 'count'
        read_count_file = f'{self.outdir_dic[sample]["mapping"]}/{sample}_raw_read_count.json'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--read_count_file {read_count_file} '

        )
        self.process_cmd(cmd, step, sample, m=2, x=1)


    def analysis(self, sample):
        step = 'analysis'
        raw_read_count_file = f'{self.outdir_dic[sample]["mapping"]}/{sample}_raw_read_count.json'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--raw_read_count_file {raw_read_count_file} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)


def main():
    multi = Multi_sweetseq(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()