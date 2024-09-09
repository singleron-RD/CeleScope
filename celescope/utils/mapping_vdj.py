from celescope.flv_trust4.__init__ import REF_DIR
from concurrent.futures import ProcessPoolExecutor
from celescope.tools import utils
import subprocess
import pysam


INDEX = {
    "TCR": ["bcrtcr", "TRA", "TRB"],
    "BCR": ["bcrtcr", "IGH", "IGK", "IGL"],
}


class Mapping_vdj:
    """
    ##Features
    - Extract reads of mapping to VDJ gene.

    For TCR library:
    Total Reads: 10,000
    Reads Mapped To Any V(D)J Genes: 9,000
    Reads Mapped To TRA: 4,000
    Reads Mapped To TRB: 5,000
    For BCR library:
    Total Reads: 20,000
    Reads Mapped To Any V(D)J Genes: 15,000
    Reads Mapped To IGH: 5,000
    Reads Mapped To IGK: 5,000
    Reads Mapped To IGL: 5,000

    ## Output
    `bcrtcr_1.fq` barcode and umi of mapping to any V(D)J genes reads.
    `bcrtcr_2.fq` sequence of mapping to any V(D)J genes reads.
    """

    def __init__(self, args):
        self.ref = args.ref
        self.fq2 = args.fq2
        self.seqtype = args.seqtype
        self.out_dir = args.out_dir
        self.out_fq1 = f"{self.out_dir}/bcrtcr_1.fq"

    def __call__(self):
        utils.check_mkdir(self.out_dir)

        results = []
        with ProcessPoolExecutor(max_workers=4) as executor:
            for result in executor.map(self.get_mapping_reads, INDEX[self.seqtype]):
                results.append(result)

        self.write_metrics()
        self.write_fq1()

    @utils.add_log
    def get_mapping_reads(self, index):
        cmd = (
            f"fastq-extractor -t 1 "
            f"-f {REF_DIR}/{self.ref}/{index}.fa "
            f"-o {self.out_dir}/{index}_2 "
            f"-u {self.fq2} "
        )
        Mapping_vdj.get_mapping_reads.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def write_metrics(self):
        metrics = []
        input_files = [self.fq2] + [
            f"{self.out_dir}/{i}_2.fq" for i in INDEX[self.seqtype]
        ]

        with ProcessPoolExecutor(max_workers=4) as executor:
            for result in executor.map(utils.get_fastx_read_number, input_files):
                metrics.append(result)

        if self.seqtype == "TCR":
            metrics = dict(
                zip(
                    [
                        "Total Reads",
                        "Reads Mapped To Any V(D)J Genes",
                        "Reads Mapped To TRA",
                        "Reads Mapped To TRB",
                    ],
                    metrics,
                )
            )
        else:
            metrics = dict(
                zip(
                    [
                        "Total Reads",
                        "Reads Mapped To Any V(D)J Genes",
                        "Reads Mapped To IGH",
                        "Reads Mapped To IGK",
                        "Reads Mapped To IGL",
                    ],
                    metrics,
                )
            )

        utils.dump_dict_to_json(metrics, f"{self.out_dir}/mapping_VDJ.json")

    @utils.add_log
    def write_fq1(self):
        out_fq1 = open(self.out_fq1, "w")
        with pysam.FastxFile(f"{self.out_dir}/bcrtcr_2.fq") as f:
            for read in f:
                bc, umi, _ = read.name.split(":")
                new_seq = bc + umi
                new_qual = "F" * len(new_seq)
                out_fq1.write(utils.fastq_line(read.name, new_seq, new_qual))
        out_fq1.close()


def mapping_vdj(args):
    runner = Mapping_vdj(args)
    runner()


def get_opts_mapping_vdj(parser, sub_program=True):
    if sub_program:
        parser.add_argument("--fq2", help="raw fastq1 file path", required=True)
        parser.add_argument(
            "--ref",
            help="reference name: hg38 or GRCm38",
            choices=["hg38", "GRCm38"],
            required=True,
        )
        parser.add_argument(
            "--seqtype", help="TCR or BCR", choices=["TCR", "BCR"], required=True
        )
        parser.add_argument("--out_dir", help="out directory", default="mapping_vdj")
    return parser
