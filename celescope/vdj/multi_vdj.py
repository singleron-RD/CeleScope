from celescope.vdj.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_vdj(Multi):

    def mapping_vdj_args(self):
        self.parser.add_argument("--type", help='TCR or BCR', required=True)
        self.parser.add_argument("--not_consensus", action='store_true', help="do not perform UMI consensus, use read instead")
        self.parser.add_argument('--species', choices=['hs', 'mmu'], help='human or mouse', default='hs')

    def read_mapping_vdj_args(self):
        self.type = self.args.type 
        self.species = self.args.species
        self.not_consensus = self.args.not_consensus
        self.not_consensus_str = Multi.arg_str(self.args.not_consensus, 'not_consensus')
    
    def count_vdj_args(self):
        self.parser.add_argument(
        '--iUMI', help='minimum number of UMI of identical receptor type and CDR3', default=1)

    def read_count_vdj_args(self):
        self.iUMI =  self.args.iUMI
    
    def custome_args(self):
        self.consensus_args()
        self.mapping_vdj_args()
        self.count_vdj_args()

    def read_custome_args(self):
        self.read_consensus_args()
        self.read_mapping_vdj_args()
        self.read_count_vdj_args()

        if self.not_consensus:
            self.steps_run = ','.join(['sample', 'barcode', 'cutadapt', 'mapping_vdj', 'count_vdj'])

    def mapping_vdj(self, sample):
        step = 'mapping_vdj'
        if self.not_consensus:
            fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        else:
            fq = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fq'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--type {self.type} '
            f'--thread {self.thread} '
            f'--species {self.species} '
            f'{self.not_consensus_str} '
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.thread)


    def count_vdj(self, sample):
        # count_vdj
        step = 'count_vdj'
        UMI_count_filter1_file = (
            f'{self.outdir_dic[sample]["mapping_vdj"]}/{sample}_UMI_count_filtered1.tsv'
        )
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--type {self.type} '
            f'--iUMI {self.iUMI} '
            f'--UMI_count_filter1_file {UMI_count_filter1_file} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=self.thread)


def main():
    multi = Multi_vdj(__ASSAY__, __STEPS__)
    multi.run()

if __name__ == '__main__':
    main()

