from celescope.vdj10X.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi
from celescope.tools.__init__ import __PATTERN_DICT__


class Multi_vdj(Multi):
    def custome_args(self):
        pass

    def read_custome_args(self):
        pass

    def convert_args(self):
        parser = self.parser
        parser.add_argument('--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version', default='auto')
        parser.add_argument('--pattern', help='')
        parser.add_argument('--whitelist', help='')
        parser.add_argument('--linker', help='')
        parser.add_argument('--lowQual', type=int, help='max phred of base as lowQual', default=10)
        parser.add_argument('--lowNum', type=int, help='max number with lowQual allowed', default=0)
        parser.add_argument('--nopolyT', action='store_true', help='output nopolyT fq')
        parser.add_argument('--noLinker', action='store_true', help='output noLinker fq')
        parser.add_argument('--probe_file', help="probe fasta file")
        parser.add_argument('--allowNoPolyT', help="allow reads without polyT", action='store_true')
        parser.add_argument('--allowNoLinker', help="allow reads without correct linker", action='store_true')
        parser.add_argument('--thread', help='number of threads', default=1)
        parser.add_argument('--method', help="read1 14bp or read2 60bp", choices=['read1', 'read2'], default='read2')
        self.parser = parser

    def read_convert_args(self):
        self.chemistry = self.args.chemistry
        self.pattern = self.args.pattern
        self.whitelist = self.args.whitelist
        self.linker = self.args.linker
        self.lowQual = self.args.lowQual
        self.lowNum = self.args.lowNum
        self.nopolyT_str = Multi.arg_str(self.args.nopolyT, 'nopolyT')
        self.noLinker_str = Multi.arg_str(self.args.noLinker, 'noLinker')
        self.probe_file = self.args.probe_file
        self.allowNoPolyT_str = Multi.arg_str(self.args.allowNoPolyT, 'allowNoPolyT')
        self.allowNoLinker_str = Multi.arg_str(self.args.allowNoLinker, 'allowNoLinker')
        self.thread = self.args.thread
        self.method = self.args.method

    def vdj_10X_args(self):
        parser = self.parser
        parser.add_argument('--species', help='species', choices=['hs','mmu'], required=True)
        parser.add_argument('--mem', help='memory (G)', default=10)
        parser.add_argument('--soft', help='cellranger version', choices=['3.0', '3.1', '6.0'], default='3.0')
        self.parser = parser

    def read_vdj_10X_args(self):
        self.mem = self.args.mem
        self.species = self.args.species
        self.soft = self.args.soft

    def convert(self, sample):
        step = 'convert'
        arr = self.fq_dict[sample]
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--chemistry {self.chemistry} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} '
            f'--pattern {self.pattern} --whitelist {self.whitelist} --linker {self.linker} '
            f'--lowQual {self.lowQual} --thread {self.thread} '
            f'--lowNum {self.lowNum} '
            f'{self.allowNoPolyT_str} '
            f'{self.allowNoLinker_str} '
            f'{self.noLinker_str} '
            f'{self.nopolyT_str} '
            f'{self.not_gzip_str} '
            f'--probe_file {self.probe_file} '
            f'--method {self.method} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


    def vdj_10X(self, sample):
        step = 'vdj_10X'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--mem {self.mem} '
            f'--thread {self.thread} '
            f'--species {self.species} '
            f'--soft {self.soft} '
        )
        self.process_cmd(cmd, step, sample, m=self.mem, x=self.thread)

    def parse_args(self):
        self.multi_opts()
        self.convert_args()
        self.vdj_10X_args()
        self.custome_args()
        self.args = self.parser.parse_args()
        # read args
        self.outdir = self.args.outdir
        self.mod = self.args.mod
        self.rm_files = self.args.rm_files
        self.steps_run = self.args.steps_run
        self.not_gzip_str = Multi.arg_str(self.args.not_gzip, 'not_gzip')
        if self.__CONDA__ == 'celescope_RD':
            self.debug_str = '--debug'
        else:
            self.debug_str = Multi.arg_str(self.args.debug, 'debug')
        self.read_convert_args()
        self.read_vdj_10X_args()
        self.read_custome_args()


def main():
    multi = Multi_vdj(__ASSAY__, __STEPS__)
    multi.run()

if __name__ == '__main__':
    main()

