from concurrent.futures import ProcessPoolExecutor
import pysam
import subprocess
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.flv_trust4.__init__ import CHAIN, INDEX, TOOLS_DIR

"""
trust goes through the following steps:
			0: start from beginning (candidate read extraction)
			1: start from assembly
			2: start from annotation
			3: start from generating the report table
"""
@utils.add_log
def extract_candidate_reads(thread, species, index_prefix, outdir, sample, fq1, fq2):
    """ Extract candiate reads for assemble.
    Args: 
        thread: number of used thread.
        species: hg38, hg19 or GRCm38.
        index_prefix: TRA, TRB or IGH, IGL, IGK.
        outdir: output directory.
        fq1: fq file including barcode and umi info.
        fq2: fq file including sequence and quality info.
    """
    cmd = (
        f'fastq-extractor -t {thread} '
        f'-f {INDEX}/{species}/{index_prefix}.fa '
        f'-o {outdir}/{sample}_{index_prefix} '
        f'--barcodeStart 0 '
        f'--barcodeEnd 23 '
        f'--umiStart 24 '
        f'--umiEnd -1 '
        f'-u {fq2} '
        f'--barcode {fq1} '
        f'--UMI {fq1} '
        )
    extract_candidate_reads.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def run_trust4(thread, fq2, fq1, IMGT_ref, fasta_ref, sample, outdir):
    """ Assemble candiate reads.
    Args: 
        thread: number of used thread.
        fq1: fq file including barcode and umi info.
        fq2: fq file including sequence and quality info.
        IMGT_ref: IMGT reference gene sequences.
        fasta_ref: V,J,C gene annotation database.
        outdir: output directory.
        sample: sample name.
    """

    cmd = (
        f'run-trust4 -t {thread} '
        f'-u {fq2} '
        f'--barcode {fq1} --barcodeRange 0 23 + '
        f'--UMI {fq1} --umiRange 24 -1 + '
        f'-f {fasta_ref} '
        f'--ref {IMGT_ref} '
        f'--outputReadAssignment -o {sample} '
        f'--od {outdir} '
    )
    run_trust4.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def get_barcode_filter_report(outdir, sample):
    """ Filter barcode report.

    If cell A's two chains CDR3s are identical to another cell B,
    and A's chain abundance is significantly lower than B's (--diffuseFrac), filter A.

    Args: 
        outdir: output directory.
        sample: sample name.
    """

    cmd = (
        f'python {TOOLS_DIR}/barcoderep-filter.py '
        f'-b {outdir}/{sample}_barcode_report.tsv > '
        f'{outdir}/{sample}_barcode_filter_report.tsv '
    )
    get_barcode_filter_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def get_filter_report(outdir, sample):
    """ Filter trust report.
    
    Filter the nonfunctional CDR3, or CDR3 sequences containing "N" in the nucleotide sequence.

    Args: 
        outdir: output directory.
        sample: sample name.
    """
    
    cmd = (f''' awk '$4!~"_" && $4!~"?"' {outdir}/{sample}_report.tsv > {outdir}/{sample}_report_filter.tsv ''')
    get_filter_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


class Assemble(Step):

    """
    ## Features

    - Assemble TCR/BCR seq data.

    ## Output
    - `03.assemble/{sample}_R1.fq` record barcode and umi info for assembly
    - `03.assemble/{sample}_R1.fq` record sequence and quality info for assembly
    - `03.assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.
    - `03.assemble/{sample}_assign.out` read assignment results.
    - `03.assemble/{sample}_assembled_reads.fa` Assembled raw reads.
    - `03.assemble/{sample}_annot.fa` Assembled annotated contig sequences.
    - `03.assemble/{sample}_report.tsv` Record assembled CDR3 types, read count and proportion of read count.
    - `03.assemble/{sample}_filter_report.tsv` Filter nonfunctional CDR3
    - `03.assemble/{sample}_barcode_report.tsv` Record chain information for each barcode.
    - `03.assemble/{sample}_barcode_filter_report.tsv` Filter low abundance cell.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        
        self.species = args.species
        self.seqtype = args.seqtype
        self.match_fq1 = args.match_fq1
        self.match_fq2 = args.match_fq2
        self.chains = self._get_chain_type(self.seqtype)
        self.matched_reads = 0


    @staticmethod
    def _get_chain_type(seqtype):
        """ Return ['TRA', 'TRB'] for TCR, ['IGH', 'IGL', 'IGK'] for BCR."""
        return CHAIN[seqtype]

    @utils.add_log
    def out_match_fastq(self):
        """
        Count matched barcodes and matched reads.
        """
        matched_cbs = set()
        with pysam.FastxFile(self.match_fq2) as fq:
            for read in fq:
                cb = read.name.split('_')[0]
                self.matched_reads += 1
                matched_cbs.add(cb)

        return matched_cbs

    @utils.add_log
    def extract_chain_reads(self):
        """ Get total and each chain type candiates read."""
        map_res = []
        map_index_prefix = self.chains
        _map_len = len(map_index_prefix)
        samples = [self.sample] * _map_len
        map_threads = [self.thread] * _map_len
        map_species = [self.species] * _map_len
        map_outdirs = [self.outdir] * _map_len
        map_fq1 = [self.match_fq1] * _map_len
        map_fq2 = [self.match_fq2] * _map_len

        with ProcessPoolExecutor(_map_len) as pool:
            for res in pool.map(extract_candidate_reads, map_threads, map_species, map_index_prefix, map_outdirs, samples, map_fq1, map_fq2):
                map_res.append(res)

    @utils.add_log
    def gen_report(self, matched_cbs):

        get_barcode_filter_report(self.outdir, self.sample)
        get_filter_report(self.outdir, self.sample)

        self.add_metric(
            name="Matched Barcodes with scRNA-seq",
            value=len(matched_cbs),
            help_info="Number of barcodes those were determined to be cells in scRNA-seq"
            )

        self.add_metric(
            name = 'Matched Reads with scRNA-seq',
            value = self.matched_reads, 
            help_info="Number of reads those were from barcodes determined to be cells in scRNA-seq"
            )

        with pysam.FastxFile(f'{self.outdir}/{self.sample}_toassemble.fq') as f:
            self.add_metric(
                name = 'Reads Mapped to Any V(D)J genes', 
                value = len(list(f)),
                total = self.matched_reads,
                help_info = "Fraction of reads that partially or wholly map to any germline V(D)J gene segment"
            )

        for _chain in self.chains:
            with pysam.FastxFile(f'{self.outdir}/{self.sample}_{_chain}.fq') as f:
                self.add_metric(
                    name = f'Reads Mapped to {_chain}', 
                    value = len(list(f)), 
                    total = self.matched_reads,
                    help_info = f"Fraction of reads that map partially or wholly to a germline {_chain} gene segment."
                )

    def run(self):
        matched_cbs = self.out_match_fastq()
        self.extract_chain_reads()

        IMGT_ref = f'{INDEX}/{self.species}/IMGT+C.fa'
        fasta_ref = f'{INDEX}/{self.species}/bcrtcr.fa'
        run_trust4(self.thread, self.match_fq2, self.match_fq1, IMGT_ref, fasta_ref, self.sample, self.outdir)
        self.gen_report(matched_cbs)



@utils.add_log
def assemble(args):
    with Assemble(args, display_title="Match and Mapping") as runner:
        runner.run()

def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--match_fq2', help='R2 reads matched with scRNA-seq.', required=True)
        parser.add_argument('--match_fq1', help='R1 reads matched with scRNA-seq.', required=True)

    parser.add_argument('--species', help='Species name and version.', choices=["hg19", "hg38", "GRCm38", "other"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--barcodeRange', help='Barcode range in fq1, INT INT CHAR.', default='0 23 +') 
    parser.add_argument('--umiRange', help='UMI range in fq1, INT INT CHAR.', default='24 -1 +')
