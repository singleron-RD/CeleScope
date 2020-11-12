from celescope.__init__ import __CONDA__
from celescope.capture_rna.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_capture_rna(Multi):
    
    def custome_args(self):
        # parser
        parser = self.parser
        parser.add_argument('--starMem', help='starMem', default=30)
        parser.add_argument('--genomeDir', help='genome index dir', required=True)
        parser.add_argument(
            '--gtf_type',
            help='Specify attribute type in GTF annotation, default=exon',
            default='exon')
        parser.add_argument('--thread', help='thread', default=6)
        parser.add_argument('--save_rds', action='store_true', help='write rds to disk')
        parser.add_argument('--type_marker_tsv', help='cell type marker tsv')
        parser.add_argument('--probe_file', help="probe fasta file")
        self.parser = parser

    def read_custome_args(self):
        self.thread = self.args.thread
        self.genomeDir = self.args.genomeDir
        self.starMem = self.args.starMem
        self.gtf_type = self.args.gtf_type
        self.save_rds = self.args.save_rds
        self.save_rds_str = Multi.arg_str(self.save_rds, 'save_rds')
        self.type_marker_tsv = self.args.type_marker_tsv
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
        self.generate_other(cmd, step, sample, m=10, x=1)


    def run_steps(self):
        for sample in self.fq_dict:
            self.sample_info(sample)
            self.barcode(sample)
            self.cutadapt(sample)
            self.STAR(sample)
            self.featureCounts(sample)
            self.count_capture_rna(sample)
            self.analysis(sample)


def main():
    multi = Multi_capture_rna(__ASSAY__, __STEPS__, __CONDA__)
    multi.col4_default = None
    multi.run()

if __name__ == '__main__':
    main()




