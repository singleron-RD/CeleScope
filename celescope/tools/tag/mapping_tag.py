import pandas as pd
import pysam

from celescope.tools import utils
from celescope.tools import parse_chemistry
from celescope.tools.step import Step, s_common


def get_opts_mapping_tag(parser, sub_program):
    parser.add_argument(
        "--fq_pattern",
        help="""R2 read pattern. The number after the letter represents the number of bases. The `fq_pattern` of CLindex is `L25C15`
`L` linker(common sequences)  
`C` tag barcode  
""",
        default="L25C15",
    )
    parser.add_argument(
        "--barcode_fasta",
        help="""Required. Tag barcode fasta file. Comma-separated if there are multiple `C` segments in fq_pattern.
The number of fasta files must equal the number of `C` segments.

It will check the mismatches between tag barcode sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < threshold. 
If no such tag exists, the read is classified as invalid.

You can find the example barcode fasta file under `celescope/data/Clindex` or `celescope/data/sweetseq`
""",
        required=True,
    )
    parser.add_argument(
        "--linker_fasta",
        help="""Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.
""",
    )
    parser.add_argument(
        "--mismatch",
        help="""Maximum number of mismatches allowed between tag barcode in R2 reads and barcode_fasta.
Comma-separated if there are multiple `C` segments (one value per segment).
If a single value is given, it applies to all `C` segments.
If not specified, it will be determined automatically for each segment based on the number of tag barcodes:
  >100000: 0
  >10000:  1
  <=10000: 2
""",
    )
    if sub_program:
        s_common(parser)
        parser.add_argument("--fq", help="R2 read fastq.", required=True)


@utils.add_log
def mapping_tag(args):
    with Mapping_tag(args, display_title="Mapping") as runner:
        runner.run()


class Mapping_tag(Step):
    """
    ## Features
    - Align R2 reads to the tag barcode fasta.

    ## Output

    - `{sample}_read_count.tsv` tab-delimited text file with 4 columns.

        `barcode` cell barcode
        `tag_name`  tag name in barcode_fasta
        `UMI`   UMI sequence
        `read_count` read count per UMI

    - `{sample}_invalid_barcode.tsv` tab-delimited text file with 2 columns.
        `tag_barcode` tag barcodes that do not match with any sequence in `--barcode_fasta`.
        `read_count` invalid tag barcode read counts
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # read args
        self.fq = args.fq
        self.fq_pattern = args.fq_pattern
        self.linker_fasta = args.linker_fasta
        self.barcode_fasta = args.barcode_fasta
        self.mismatch_arg = args.mismatch

        # process pattern
        self.pattern_dict = parse_chemistry.parse_pattern(self.fq_pattern)
        self.c_slices = self.pattern_dict["C"]
        self.n_c = len(self.c_slices)

        # parse barcode_fasta
        fasta_list = [f.strip() for f in self.barcode_fasta.split(",")]
        if len(fasta_list) != self.n_c:
            raise ValueError(
                f"Number of barcode_fasta files ({len(fasta_list)}) "
                f"must equal number of `C` segments in fq_pattern ({self.n_c})."
            )

        # read barcode fastas
        self.barcode_dict_list = []
        self.barcode_length_list = []
        for i, fasta in enumerate(fasta_list):
            bc_dict, bc_len = utils.read_fasta(fasta, equal=True)
            expected_len = self.c_slices[i].stop - self.c_slices[i].start
            if bc_len != expected_len:
                raise ValueError(
                    f"Length of tag barcode in fasta[{i}] ({bc_len}) != "
                    f"length of C[{i}] in fq_pattern ({expected_len})."
                )
            self.barcode_dict_list.append(bc_dict)
            self.barcode_length_list.append(bc_len)

        # parse mismatch
        self.n_mismatch_list = self._parse_mismatch()

        # linker
        if self.linker_fasta and self.linker_fasta != "None":
            self.linker_dict, self.linker_length = utils.read_fasta(
                self.linker_fasta, equal=True
            )
            len_L = sum(x.stop - x.start for x in self.pattern_dict["L"])
            if len_L != self.linker_length:
                raise ValueError(f"""The length of linker in fq_pattern({len_L}) != 
                    length of linker in linker_fasta({self.linker_length})""")
        else:
            self.linker_dict, self.linker_length = {}, 0

        # mismatch dicts per segment
        self.mismatch_dict_list = self._get_tag_barcode_mismatch_dicts()

        # variables
        self.total_reads = 0
        self.reads_unmapped_too_short = 0
        self.reads_unmapped_invalid_linker = 0
        self.reads_unmapped_invalid_barcode = 0
        self.reads_mapped = 0
        self.res_dic = utils.nested_defaultdict()
        self.res_sum_dic = utils.nested_defaultdict(dim=2)
        self.match_barcode = []
        self.invalid_barcode_dict = utils.nested_defaultdict(dim=1)

        # out files
        self.read_count_file = f"{self.outdir}/{self.sample}_read_count.tsv"
        self.invalid_barcode_file = f"{self.outdir}/{self.sample}_invalid_barcode.tsv"

    def _parse_mismatch(self):
        if self.mismatch_arg is not None:
            mismatch_strs = [m.strip() for m in self.mismatch_arg.split(",")]
            mismatch_vals = [int(m) for m in mismatch_strs]
            if len(mismatch_vals) == 1 and self.n_c > 1:
                mismatch_vals = mismatch_vals * self.n_c
            if len(mismatch_vals) != self.n_c:
                raise ValueError(
                    f"Number of mismatch values ({len(mismatch_vals)}) "
                    f"must equal number of `C` segments ({self.n_c})."
                )
        else:
            mismatch_vals = []
            for bc_dict in self.barcode_dict_list:
                n_barcodes = len(bc_dict)
                if n_barcodes > 100000:
                    mismatch_vals.append(0)
                elif n_barcodes > 10000:
                    mismatch_vals.append(1)
                else:
                    mismatch_vals.append(2)

        for i in range(self.n_c):
            if mismatch_vals[i] > self.barcode_length_list[i]:
                mismatch_vals[i] = self.barcode_length_list[i]

        return mismatch_vals

    @utils.add_log
    def _get_tag_barcode_mismatch_dicts(self):
        mismatch_dict_list = []
        for i in range(self.n_c):
            mismatch_dict = {}
            n_mismatch = self.n_mismatch_list[i]
            for seq_id, seq in self.barcode_dict_list[i].items():
                for mismatch_seq in parse_chemistry.create_mismatch_seqs(
                    seq, max_mismatch=n_mismatch
                ):
                    mismatch_dict[mismatch_seq] = seq_id
            mismatch_dict_list.append(mismatch_dict)
        return mismatch_dict_list

    def check_barcode_with_mismatch(self, barcode, seq_barcode_list, umi):
        """
        Args:
            barcode: cell barcode
            seq_barcode_list: list of tag barcode sequences, one per C segment
            umi: UMI sequence
        """
        matched_names = []
        for i, seq_barcode in enumerate(seq_barcode_list):
            if seq_barcode not in self.mismatch_dict_list[i]:
                self.reads_unmapped_invalid_barcode += 1
                self.invalid_barcode_dict[seq_barcode] += 1
                return
            matched_names.append(self.mismatch_dict_list[i][seq_barcode])

        tag_name = "_".join(matched_names)
        self.res_dic[barcode][tag_name][umi] += 1
        self.reads_mapped += 1

    def process_read(self):
        with pysam.FastxFile(self.fq) as infile:
            for record in infile:
                self.total_reads += 1
                attr = str(record.name).strip("@").split(":")
                barcode = str(attr[0])
                umi = str(attr[1])
                seq = record.sequence

                if self.linker_length != 0:
                    seq_linker = "".join(seq[x] for x in self.pattern_dict["L"])

                seq_barcode_list = [seq[s] for s in self.c_slices]

                # check linker
                if self.linker_length != 0:
                    valid_linker = False
                    for linker_name in self.linker_dict:
                        if utils.hamming_correct(
                            self.linker_dict[linker_name], seq_linker
                        ):
                            valid_linker = True
                            break
                else:
                    valid_linker = True

                if not valid_linker:
                    self.reads_unmapped_invalid_linker += 1
                    continue

                # check barcode per segment
                self.check_barcode_with_mismatch(barcode, seq_barcode_list, umi)

    def write_files(self):
        # write dic to pandas df
        rows = []
        for barcode in self.res_dic:
            for tag_name in self.res_dic[barcode]:
                for umi in self.res_dic[barcode][tag_name]:
                    rows.append(
                        [barcode, tag_name, umi, self.res_dic[barcode][tag_name][umi]]
                    )
        df_read_count = pd.DataFrame(rows)
        df_read_count.rename(
            columns={0: "barcode", 1: "tag_name", 2: "UMI", 3: "read_count"},
            inplace=True,
        )
        df_read_count.to_csv(self.read_count_file, sep="\t", index=False)

        # write invalid seq_barcode to file
        with open(self.invalid_barcode_file, "w") as f:
            f.write("tag_barcode\tread_count\n")
            for seq_barcode, count in sorted(
                self.invalid_barcode_dict.items(), key=lambda x: x[1], reverse=True
            ):
                f.write(f"{seq_barcode}\t{count}\n")

    def add_metrics(self):
        # add metrics
        self.add_metric(
            name="Reads Mapped",
            value=self.reads_mapped,
            total=self.total_reads,
            help_info="R2 reads that successfully mapped to linker and tag-barcode",
        )
        self.add_metric(
            name="Reads Unmapped Invalid Linker",
            value=self.reads_unmapped_invalid_linker,
            total=self.total_reads,
            help_info="Unmapped R2 reads because of too many mismatches in linker sequence",
        )
        self.add_metric(
            name="Reads Unmapped Invalid Barcode",
            value=self.reads_unmapped_invalid_barcode,
            total=self.total_reads,
            help_info="Unmapped R2 reads because of too many mismatches in tag-barcode sequence",
        )

    @utils.add_log
    def run(self):
        self.process_read()
        self.write_files()
        self.add_metrics()
