import os
import pathlib
from collections import defaultdict
from itertools import groupby
import subprocess
import unittest
import shutil

import pysam

from celescope.rna.mkref import Mkref_rna
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.tools import reference
from celescope.tools.__init__ import TAG_BAM_SUFFIX


GTF_TYPES = ["exon", "gene"]


def correct_umi(umi_dict, percent=0.1):
    """
    Correct umi_dict in place.
    Args:
        umi_dict: {umi_seq: umi_count}
        percent: if hamming_distance(low_seq, high_seq) == 1 and
            low_count / high_count < percent, merge low to high.
        return_dict: if set, return correct_dict = {low_seq:high_seq}
    Returns:
        n_corrected_umi: int
        n_corrected_read: int
    """
    n_corrected_umi = 0
    n_corrected_read = 0
    dic = {}

    # sort by value(UMI count) first, then key(UMI sequence)
    umi_arr = sorted(umi_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    while True:
        # break when only highest in umi_arr
        if len(umi_arr) == 1:
            break
        umi_low = umi_arr.pop()
        low_seq = umi_low[0]
        low_count = umi_low[1]

        for umi_kv in umi_arr:
            high_seq = umi_kv[0]
            high_count = umi_kv[1]
            if float(low_count / high_count) > percent:
                break
            if utils.hamming_distance(low_seq, high_seq) == 1:
                n_low = umi_dict[low_seq]
                n_corrected_umi += 1
                n_corrected_read += n_low
                # merge
                umi_dict[high_seq] += n_low
                dic[low_seq] = high_seq
                del umi_dict[low_seq]
                break

    return n_corrected_umi, n_corrected_read, dic


def discard_read(gene_umi_dict):
    """
    If two or more groups of reads have the same barcode and UMI, but different gene annotations, the gene annotation with the most supporting reads is kept for UMI counting, and the other read groups are discarded. In case of a tie for maximal read support, all read groups are discarded, as the gene cannot be confidently assigned.

    Returns:
        discarded_umi: set. umi with tie read count
        umi_gene_dict: {umi_seq: {gene_id: read_count}}
    """

    discard_umi = set()
    umi_gene_dict = defaultdict(lambda: defaultdict(int))
    for gene_id in gene_umi_dict:
        for umi in gene_umi_dict[gene_id]:
            umi_gene_dict[umi][gene_id] += gene_umi_dict[gene_id][umi]

    for umi in umi_gene_dict:
        max_read_count = max(umi_gene_dict[umi].values())
        gene_id_max = [
            gene_id
            for gene_id, read_count in umi_gene_dict[umi].items()
            if read_count == max_read_count
        ]

        if len(gene_id_max) > 1:
            discard_umi.add(umi)
        else:
            gene_id = gene_id_max[0]
            umi_gene_dict[umi] = {gene_id: umi_gene_dict[umi][gene_id]}

    return discard_umi, umi_gene_dict


class FeatureCounts(Step):
    """
    ## Features
    - Assigning uniquely mapped reads to genomic features with FeatureCounts.
    ## Output
    - `{sample}` Numbers of reads assigned to features (or meta-features).
    - `{sample}_summary` Stat info for the overall summrization results, including number of
    successfully assigned reads and number of reads that failed to be assigned due to
    various reasons (these reasons are included in the stat info).
    - `{sample}_aligned_sortedByCoord_addTag.bam` featureCounts output BAM,
    sorted by coordinates
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # set
        self.gtf = Mkref_rna.get_config(self.args.genomeDir)["files"]["gtf"]
        gp = reference.GtfParser(self.gtf)
        self.id_name = gp.get_id_name()
        self.intron_dict = {}

        # stats
        self.exon = self.intron = self.intergenic = self.ambiguity = 0

        # temp file
        self.tmp_dir = f"{self.outdir}/tmp/"
        self.add_tag_bam = f"{self.out_prefix}_addTag.bam"
        input_basename = os.path.basename(self.args.input)
        self.exon_bam = f"{self.tmp_dir}/exon/{input_basename}.featureCounts.bam"
        self.intron_bam = f"{self.tmp_dir}/intron/{input_basename}.featureCounts.bam"

        # out
        self.count_detail_file = f"{self.out_prefix}_count_detail.txt"
        self.nameSorted_bam = f"{self.out_prefix}_nameSorted.bam"
        self.out_bam = f"{self.out_prefix}_{TAG_BAM_SUFFIX}"

    def add_tag(self, seg, id_name):
        """
        Add intron reads and tag

        Args:
            seg: pysam bam segment
            id_name: {gene_id: gene_name}

        Returns:
            seg with tag added

        Tags:
            CB: cell barcode
            UB: error-corrected UMI
            UR: original UMI
            GN: gene name
            GX: gene_id

        """
        attr = seg.query_name.split(":")
        barcode = attr[0]
        ur = ub = attr[1]

        # assign to some gene
        xs = seg.get_tag("XS")
        if xs == "Assigned":
            gene_id = seg.get_tag("XT")
            gene_name = id_name[gene_id]
            seg.set_tag(tag="GN", value=gene_name, value_type="Z")
            seg.set_tag(tag="GX", value=gene_id, value_type="Z")
            seg.set_tag(tag="RE", value="E", value_type="Z")
            self.exon += 1
        else:
            if self.intron_dict and seg.query_name in self.intron_dict:
                gene_id = self.intron_dict[seg.query_name]
                gene_name = id_name[gene_id]
                seg.set_tag(tag="GN", value=gene_name, value_type="Z")
                seg.set_tag(tag="GX", value=gene_id, value_type="Z")
                seg.set_tag(tag="RE", value="N", value_type="Z")
                seg.set_tag(tag="XT", value=gene_id, value_type="Z")
                self.intron += 1
            elif xs == "Unassigned_NoFeatures":
                seg.set_tag(tag="RE", value="I", value_type="Z")
                self.intergenic += 1
            elif xs == "Unassigned_Ambiguity":
                seg.set_tag(tag="RE", value="A", value_type="Z")
                self.ambiguity += 1

        seg.set_tag(tag="CB", value=barcode, value_type="Z")
        seg.set_tag(tag="UB", value=ub, value_type="Z")
        seg.set_tag(tag="UR", value=ur, value_type="Z")

        return seg

    @utils.add_log
    def get_intron_dict(self):
        with pysam.AlignmentFile(self.intron_bam, "rb") as in_bam:
            for seg in in_bam:
                if seg.has_tag("XT"):
                    self.intron_dict[seg.query_name] = seg.get_tag("XT")

    @utils.add_log
    def run_featureCounts(self, outdir, gtf_type):
        """
        allow multimapping with -M; but each multi-mapped reads only have one alignment because of --outSAMmultNmax 1
        """
        cmd = (
            "featureCounts "
            f"-s 1 "
            f"--largestOverlap "
            f"-M "
            f"-a {self.gtf} "
            f"-o {outdir}/{self.sample} "
            "-R BAM "
            f"-T {self.args.thread} "
            f"-t {gtf_type} "
            f"{self.args.input} "
            "2>&1 "
        )
        if self.args.featureCounts_param:
            cmd += " " + self.args.featureCounts_param
        self.run_featureCounts.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run_exon_intron(self):
        tmp_dir = f"{self.outdir}/tmp/"
        for gtf_type in ["exon", "intron"]:
            outdir = f"{tmp_dir}/{gtf_type}"
            pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
            self.run_featureCounts(outdir, gtf_type)

    def remove_temp_file(self):
        shutil.rmtree(self.tmp_dir)
        os.remove(self.add_tag_bam)

    def run(self):
        self.run_exon_intron()
        self.get_intron_dict()
        utils.sort_bam(
            self.exon_bam, self.nameSorted_bam, threads=self.thread, by="name"
        )
        self.get_count_detail_add_tag()
        self.add_metrics()
        utils.sort_bam(input_bam=self.add_tag_bam, output_bam=self.out_bam)
        self.remove_temp_file()

    @utils.add_log
    def add_metrics(self):
        total = self.exon + self.intron + self.intergenic + self.ambiguity

        self.add_metric(
            name="Feature Type",
            value=self.args.gtf_type.capitalize(),
            help_info="Specified by `--gtf_type`. For snRNA-seq, you need to add `--gtf_type gene` to include reads mapped to intronic regions. Staring from CeleScope v1.12.0, the default value of gtf_type is changed from `exon` to `gene`.",
        )
        self.add_metric(
            name="Reads Assigned To Exonic Regions",
            value=self.exon,
            total=total,
            help_info="Reads that can be successfully assigned to exonic regions",
        )
        self.add_metric(
            name="Reads Assigned To Intronic Regions",
            value=self.intron,
            total=total,
            help_info="Reads that can be successfully assigned to intronic regions",
        )
        self.add_metric(
            name="Reads Assigned To Intergenic Regions",
            value=self.intergenic,
            total=total,
            help_info="Reads that can be successfully assigned to intergenic regions",
        )
        self.add_metric(
            name="Reads Unassigned Ambiguity",
            value=self.ambiguity,
            total=total,
            help_info="Alignments that overlap two or more features",
        )

    @utils.add_log
    def get_count_detail_add_tag(self):
        """
        bam to detail table
        must be used on name_sorted bam
        Output file:
            - count_detail_file
            - bam with tag(remain name sorted)
        """
        save = pysam.set_verbosity(0)
        inputFile = pysam.AlignmentFile(self.nameSorted_bam, "rb")
        outputFile = pysam.AlignmentFile(
            self.add_tag_bam, "wb", header=inputFile.header
        )
        pysam.set_verbosity(save)

        with open(self.count_detail_file, "wt") as fh1:
            fh1.write(
                "\t".join(
                    ["Barcode", "geneID", "UMI", "read", "unique", "PCR_duplicate"]
                )
                + "\n"
            )

            def keyfunc(x):
                return x.query_name.split(":", 1)[0]

            for _, g in groupby(inputFile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                gene_umi_pos = utils.nested_defaultdict(dim=3, valType=int)
                for seg in g:
                    seg = self.add_tag(seg, self.id_name)
                    outputFile.write(seg)
                    barcode, umi = seg.get_tag("CB"), seg.get_tag("UB")
                    if not seg.has_tag("GX"):
                        continue
                    gene_id = seg.get_tag("GX")
                    gene_umi_dict[gene_id][umi] += 1
                    gene_umi_pos[gene_id][umi][seg.reference_start] += 1

                # output
                for gene_id in gene_umi_dict:
                    n_umi = len(gene_umi_dict[gene_id])
                    n_read = 0
                    unique = dup = 0
                    for umi in gene_umi_dict[gene_id]:
                        read_count = gene_umi_dict[gene_id][umi]
                        n_read += read_count
                        if read_count == 1:
                            # unique
                            unique += 1
                        else:
                            # only add postion duplicate read number
                            for pos in gene_umi_pos[gene_id][umi]:
                                if gene_umi_pos[gene_id][umi][pos] > 1:
                                    dup += gene_umi_pos[gene_id][umi][pos]
                    fh1.write(
                        f"{barcode}\t{gene_id}\t{n_umi}\t{n_read}\t{unique}\t{dup}\n"
                    )

        inputFile.close()
        outputFile.close()


@utils.add_log
def featureCounts(args):
    with FeatureCounts(args) as runner:
        runner.run()


def get_opts_featureCounts(parser, sub_program):
    parser.add_argument(
        "--gtf_type",
        help="Specify feature type in GTF annotation",
        default="gene",
        choices=["exon", "gene"],
    )
    parser.add_argument("--genomeDir", help=HELP_DICT["genomeDir"])
    parser.add_argument(
        "--featureCounts_param", help=HELP_DICT["additional_param"], default=""
    )
    # parser.add_argument('--correct_UMI', help='perform UMI correction.')

    if sub_program:
        parser.add_argument("--input", help="Required. BAM file path.", required=True)
        parser = s_common(parser)
    return parser


class featureCounts_test(unittest.TestCase):
    def test_correct_umi(self):
        dic = {
            "apple1": 2,
            "apple2": 30,
            "bears1": 5,
            "bears2": 10,
            "bears3": 100,
            "ccccc1": 20,
            "ccccc2": 199,
        }
        n_corrected_umi, n_corrected_read, _ = correct_umi(dic)
        dic_after_correct = {
            "ccccc1": 20,
            "apple2": 32,
            "bears3": 115,
            "ccccc2": 199,
        }
        self.assertEqual(dic, dic_after_correct)
        self.assertEqual(n_corrected_umi, 3)
        self.assertEqual(n_corrected_read, 2 + 5 + 10)


if __name__ == "__main__":
    unittest.main()
