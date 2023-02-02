from multiprocessing import Pool
import subprocess
import math

from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.flv_trust4.__init__ import CHAIN, REF_DIR

CANDIDATE_FQ_SUFFIX = 'bcrtcr.fq'


class Mapping(Step):
    """
    ## Features

    - Extract candidate reads to assemble.

    ## Output
    - `02.mapping/{sample}_bcrtcr.fq` All candidate reads(mapped to any V(D)J genes) sequence.
    - `02.mapping/{sample}_bcrtcr_bc.fa` All candidate reads(mapped to any V(D)J genes) barcode.
    - `02.mapping/{sample}_bcrtcr_umi.fa` All candidate reads(mapped to any V(D)J genes) umi.
    - `02.mapping/{sample}_{chain}.fq` Candidate reads(mapped to {chain}) sequence. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.
    - `02.mapping/{sample}_{chain}_bc.fa` Candidate reads(mapped to {chain}) barcode. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.
    - `02.mapping/{sample}_{chain}_umi.fa` Candidate reads(mapped to {chain}) umi. For BCR: IGH, IGK, IGL, For TCR: TRA, TRB.
    """
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self._ref = args.ref
        self._barcodeRange = args.barcodeRange
        self._umiRange = args.umiRange
        self._match_fq1 = args.match_fq1
        self._match_fq2 = args.match_fq2
        self._chains = CHAIN[args.seqtype]

        n_extract = len(self._chains) + 1
        self._single_thread = math.ceil(self.thread / n_extract)
        self._matched_reads = self.get_slot_key(
            slot='metrics',
            step_name='barcode',
            key='Valid Matched Reads',
        )

    
    @utils.add_log
    def extract_chain_reads(self):
        """
        bcrtcr.fq + {chains}.fq
        """
        cb_range = self._barcodeRange.split(' ')
        umi_range = self._umiRange.split(' ')

        map_index_prefix = ['bcrtcr'] + self._chains
        n_map = len(map_index_prefix)
        samples = [self.sample] * n_map
        map_ref = [self._ref] * n_map
        map_outdirs = [self.outdir] * n_map
        map_fq1 = [self._match_fq1] * n_map
        map_fq2 = [self._match_fq2] * n_map
        map_cb_range = [cb_range] * n_map
        map_umi_range = [umi_range] * n_map
        map_n_thread = [self._single_thread] * n_map

        with Pool(n_map) as pool:
            pool.starmap(Mapping.extract_candidate_reads, 
            zip(map_ref, map_index_prefix, map_outdirs, samples, map_fq1, map_fq2, map_cb_range, map_umi_range, map_n_thread))  

    @staticmethod
    @utils.add_log
    def extract_candidate_reads(ref, index_prefix, outdir, sample, fq1, fq2, barcodeRange, umiRange, single_thread):
        """
        helper function for extract reads map to index_prefix
        index_prefix can be 'bcrtcr' + chains
        """
        cmd = (
            f'fastq-extractor -t {single_thread} '
            f'-f {REF_DIR}/{ref}/{index_prefix}.fa '
            f'-o {outdir}/{sample}_{index_prefix} '
            f'--barcodeStart {barcodeRange[0]} '
            f'--barcodeEnd {barcodeRange[1]} '
            f'--umiStart {umiRange[0]} '
            f'--umiEnd {umiRange[1]} '
            f'-u {fq2} '
            f'--barcode {fq1} '
            f'--UMI {fq1} '
            '2>&1 '
            )
        Mapping.extract_candidate_reads.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)  

    def add_metrics(self):
        n_bcrtcr = utils.get_fastx_read_number(f'{self.out_prefix}_{CANDIDATE_FQ_SUFFIX}')
        self.add_metric(
            name = 'Reads Mapped to Any V(D)J genes', 
            value = n_bcrtcr,
            total = self._matched_reads,
            help_info = "Fraction of reads that partially or wholly map to any germline V(D)J gene segment"
        )

        for _chain in self._chains:
            n_chain = utils.get_fastx_read_number(f'{self.out_prefix}_{_chain}.fq')
            self.add_metric(
                name = f'Reads Mapped to {_chain}', 
                value = n_chain, 
                total = self._matched_reads,
                help_info = f"Fraction of reads that map partially or wholly to a germline {_chain} gene segment."
            )

    def run(self):
        self.extract_chain_reads()
        self.add_metrics()


@utils.add_log
def mapping(args):
    with Mapping(args, display_title="Mapping") as runner:
        runner.run()


def get_opts_mapping(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--match_fq2', help='R2 reads matched with scRNA-seq.', required=True)
        parser.add_argument('--match_fq1', help='R1 reads matched with scRNA-seq.', required=True)

    parser.add_argument('--ref', help='reference name', choices=["hg19", "hg38", "GRCm38", "other"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--barcodeRange', help='Barcode range in fq1, INT INT CHAR.', default='0 23 +') 
    parser.add_argument('--umiRange', help='UMI range in fq1, INT INT CHAR.', default='24 -1 +')