import subprocess
from celescope.tools.step import Step, s_common


class PathSeq(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.unmap_bam = f"{self.out_prefix}_unmapped.bam"
        self.addRG_bam = f"{self.out_prefix}_addRG.bam"
        self.pathseq_bam = f"{self.out_prefix}_pathseq.bam"
        self.pathseq_score_txt = f"{self.out_prefix}_pathseq_score.txt"
        self.filter_metrics_txt = f"{self.out_prefix}_filter_metrics.txt"

    def get_unmap_bam(self):
        cmd = (
            f"samtools view -b -h -f 4 "
            f"{self.args.input_bam} "
            f"> {self.unmap_bam} "
        )
        subprocess.check_call(cmd, shell=True)

    def picard_addRG(self):
        cmd = (
            f"picard AddOrReplaceReadGroups "
            f"INPUT={self.unmap_bam} "
            f"OUTPUT={self.addRG_bam} "
            f"RGID={self.sample}_rgid "
            f"RGLB={self.sample}_rglb "
            f"RGPL=illumina "
            f"RGPU={self.sample}_rgpu "
            f"RGSM={self.sample}_rgsm "
        )
        subprocess.check_call(cmd, shell=True)

    def run_pathseq(self):
        cmd = (
            f"gatk PathSeqPipelineSpark "
            f"--input {self.addRG_bam} "
            f"--filter-bwa-image {self.args.filter_bwa_image} "
            f"--kmer-file {self.args.kmer_file} "
            f"--microbe-bwa-image {self.args.microbe_bwa_image} "
            f"--microbe-dict {self.args.microbe_dict} "
            f"--taxonomy-file {self.args.microbe_taxonomy_file} "
            f"--java-options -Xmx100g "
            f"--min-score-identity {self.args.min_score_identity} "
            f"--min-clipped-read-length {self.args.min_clipped_read_length} "
            f"--is-host-aligned true "
            f"--filter-duplicates false "
            f"--output {self.pathseq_bam} "
            f"--scores-output {self.pathseq_score_txt} "
            f"--filter-metrics {self.filter_metrics_txt} "
        )
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.get_unmap_bam()
        self.picard_addRG()
        self.run_pathseq()


def pathseq(args):
    with PathSeq(args) as runner:
        runner.run()


def get_opts_pathseq(parser, sub_program):
    parser.add_argument("--filter_bwa_image", help="filter_bwa_image", required=True)
    parser.add_argument("--kmer_file", help="kmer_file", required=True)
    parser.add_argument("--microbe_bwa_image", help="microbe_bwa_image", required=True)
    parser.add_argument("--microbe_dict", help="microbe_dict", required=True)
    parser.add_argument(
        "--microbe_taxonomy_file", help="microbe_taxonomy_file", required=True
    )
    parser.add_argument("--min_score_identity", help="min_score_identity", default=0.7)
    parser.add_argument(
        "--min_clipped_read_length", help="min_clipped_read_length", default=60
    )
    if sub_program:
        parser.add_argument("--input_bam", help="input bam file", required=True)
        parser = s_common(parser)
    return parser
