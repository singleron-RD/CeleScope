from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import pysam
import subprocess
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.fl_vdj_TRUST4.__init__ import CHAIN
from celescope.fl_vdj_TRUST4.__init__ import INDEX, TOOLS_DIR

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
    
    cmd = (f''' awk '$4!~"_" && $4!~"?"' {outdir}/{sample}_report.tsv > {outdir}/{sample}_filter_report.tsv ''')
    get_filter_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def get_full_len_assembly(outdir, sample):
    """ Get full length assembled contig.

    get full length result from annotation of the consensus assembly in fasta format.

    Args: 
        outdir: output directory.
        sample: sample name.
    """

    cmd = (
        f'perl {TOOLS_DIR}/GetFullLengthAssembly.pl '
        f'{outdir}/{sample}_annot.fa > '
        f'{outdir}/{sample}_full_len.fa '
    )
    get_full_len_assembly.logger.info(cmd)
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
    - `03.assemble/{sample}_full_len.fa` Assembled full length contig sequences.
    - `03.assemble/{sample}_report.tsv` Record assembled CDR3 types, read count and proportion of read count.
    - `03.assemble/{sample}_filter_report.tsv` Filter nonfunctional CDR3
    - `03.assemble/{sample}_barcode_report.tsv` Record chain information for each barcode.
    - `03.assemble/{sample}_barcode_filter_report.tsv` Filter low abundance cell.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        
        self.cutadapted_fq = args.cutadapted_fq
        self.species = args.species
        self.seqtype = args.seqtype
        self.match_dir = args.match_dir
        self.match_previous_assemble = args.match_previous_assemble
        self.fq1, self.fq2 = f'{self.outdir}/{self.sample}_R1.fq', f'{self.outdir}/{self.sample}_R2.fq'
        
        self.chains = self._get_chain_type(self.seqtype)
        self.match_barcodes = self._get_match_barcode(self.match_dir)

    @staticmethod
    def reversed_compl(seq):
        """ Return reverse complementation sequence."""
        return str(Seq(seq).reverse_complement())

    @staticmethod
    def _get_chain_type(seqtype):
        """ Return ['TRA', 'TRB'] for TCR, ['IGH', 'IGL', 'IGK'] for BCR."""
        return CHAIN[seqtype]
    
    @staticmethod
    def _get_match_barcode(match_dir):
        """ Return matched barcodes set for match scRNA-seq dir."""
        match_barcodes, _ = utils.get_barcode_from_match_dir(match_dir) 
        return match_barcodes

    @utils.add_log
    def gen_fq(self):
        """ Generate fq1 and fq2 file.
        
        fq1 file including barcode and umi info.
        fq2 file including sequence and quality info.

        Returns:
            matched_cbs: set of matched barcodes with scRNA-seq.
            matched_read_count: read count for matched barcode with scRNA-seq.
            read_count: read count for all barcode.
        """

        out_fq1 = open(self.fq1, 'w')
        out_fq2 = open(self.fq2, 'w')
        read_dict = defaultdict(list)

        total_read_count, matched_read_count = 0, 0
        with pysam.FastxFile(self.cutadapted_fq) as fq:
            for read in fq:
                attr = read.name.split('_')
                bc = attr[0]
                read_dict[bc].append(read)
            barcodes = set(read_dict.keys())

            rna_cbs = [self.reversed_compl(cb) for cb in self.match_barcodes]
            matched_cbs = set(read_dict.keys()).intersection(set(rna_cbs))

            if self.match_previous_assemble:
                barcodes = matched_cbs

            assert len(barcodes) != 0

            for cb in barcodes:
                for read in read_dict[cb]:
                    umi = read.name.split('_')[1]
                    qual = 'F' * len(cb + umi)
                    seq1 = f'@{read.name}\n{cb}{umi}\n+\n{qual}\n'
                    out_fq1.write(seq1)
                    out_fq2.write(str(read)+'\n')
                    total_read_count += 1
                    if cb in matched_cbs:
                        matched_read_count += 1
            out_fq1.close()
            out_fq2.close()

            read_count = total_read_count
            if self.match_previous_assemble:
                read_count = matched_read_count
            
        return matched_cbs, matched_read_count, read_count

    @utils.add_log
    def get_VDJmapping_reads(self):
        """ Get total and each chain type candiates read."""
        map_res = []
        map_index_prefix = self.chains
        _map_len = len(map_index_prefix)
        samples = [self.sample] * _map_len
        map_threads = [self.thread] * _map_len
        map_species = [self.species] * _map_len
        map_outdirs = [self.outdir] * _map_len
        fq1 = [self.fq1] * _map_len
        fq2 = [self.fq2] * _map_len

        with ProcessPoolExecutor(_map_len) as pool:
            for res in pool.map(extract_candidate_reads, map_threads, map_species, map_index_prefix, map_outdirs, samples, fq1, fq2):
                map_res.append(res)

    @utils.add_log
    def gen_report(self, matched_cbs, matched_read_count, read_count):
        get_barcode_filter_report(self.outdir, self.sample)
        get_filter_report(self.outdir, self.sample)
        self.get_VDJmapping_reads()

        self.add_metric(
            name="Matched Barcodes with scRNA-seq",
            value=len(matched_cbs),
            help_info="Number of barcodes those were determined to be cells in scRNA-seq"
            )

        self.add_metric(
            name = 'Matched Reads with scRNA-seq',
            value = matched_read_count, 
            help_info="Number of reads those were from barcodes determined to be cells in scRNA-seq"
            )

        with pysam.FastxFile(f'{self.outdir}/{self.sample}_toassemble.fq') as f:
            self.add_metric(
                name = 'Reads Mapped to Any V(D)J genes', 
                value = len(list(f)),
                total = read_count,
                help_info = "Fraction of reads that partially or wholly map to any germline V(D)J gene segment"
            )

        for _chain in self.chains:
            with pysam.FastxFile(f'{self.outdir}/{self.sample}_{_chain}.fq') as f:
                self.add_metric(
                    name = f'Reads Mapped to {_chain}', 
                    value = len(list(f)), 
                    total = read_count,
                    help_info = f"Fraction of reads that map partially or wholly to a germline {_chain} gene segment."
                )

    def run(self):
        matched_cbs, matched_read_count, read_count = self.gen_fq()

        IMGT_ref = f'{INDEX}/{self.species}/IMGT+C.fa'
        fasta_ref = f'{INDEX}/{self.species}/bcrtcr.fa'
        run_trust4(self.thread, self.fq2, self.fq1, IMGT_ref, fasta_ref, self.sample, self.outdir)

        self.gen_report(matched_cbs, matched_read_count, read_count)
        get_full_len_assembly(self.outdir, self.sample)


@utils.add_log
def assemble(args):
    with Assemble(args, display_title="Match and Mapping") as runner:
        runner.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--cutadapted_fq', help='cutadapted_fq', required=True)
        parser.add_argument('--match_dir', help='Match scRNA-seq directory.', required=True)

    parser.add_argument('--species', help='Species name and version.', choices=["hg19", "hg38", "GRCm38"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--match_previous_assemble', help='whether match reads with sc-RNA before assemble', action='store_true')