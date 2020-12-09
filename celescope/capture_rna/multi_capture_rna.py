from celescope.__init__ import __CONDA__
from celescope.capture_rna.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_capture_rna(Multi):
    
    def custome_args(self):
        # parser
        self.STAR_args()
        self.analysis_args()
        self.parser.add_argument('--probe_file', help="probe fasta file")

    def read_custome_args(self):
        self.read_STAR_args()
        self.read_analysis_args()
        self.probe_file = self.args.probe_file

    def count_capture_rna(self, sample):
        step = 'count_capture_rna'
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_name_sorted.bam'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step } '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--bam {bam} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--genomeDir {self.genomeDir} '
            f'--cells auto '
        )
        self.process_cmd(cmd, step, sample, m=10, x=1)
    
    def analysis(self, sample):
        step = 'analysis'
        matrix_file = f'{self.outdir_dic[sample]["count_capture_rna"]}/{sample}_matrix.tsv.gz'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step } '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--matrix_file {matrix_file} '
            f'{self.save_rds_str} '
            f'--type_marker_tsv {self.type_marker_tsv} '
        )
        self.process_cmd(cmd, step, sample, m=10, x=1)


def main():
    multi = Multi_capture_rna(__ASSAY__, __STEPS__, __CONDA__)
    multi.col4_default = None
    multi.run()

if __name__ == '__main__':
    main()




