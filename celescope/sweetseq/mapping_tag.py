from celescope.tools.tag.mapping_tag import Mapping_tag as Mt, mapping_tag as mt
from celescope.tools.step import s_common


# L23C15, not L25C15 like Clindex
def get_opts_mapping_tag(parser, sub_program):
    parser.add_argument(
        "--fq_pattern",
        help="""R2 read pattern. The number after the letter represents the number of bases. CLindex is `L25C15` and sweetseq is `L23C15`
`L` linker(common sequences)  
`C` tag barcode  
""",
        default="L23C15",
    )
    parser.add_argument(
        "--barcode_fasta",
        help="""Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < threshold. 
If no such tag exists, the read is classified as invalid.

You can find the barcode fasta file under `celescope/data/sweetseq`
""",
        required=True,
    )
    parser.add_argument(
        "--linker_fasta",
        help="""Optional. If provided, it will check the mismatches between linker sequence in R2 reads 
with all linker sequence in linker_fasta. If no mismatch < len(linker) / 10 + 1, the read is classified as invalid.
""",
    )
    if sub_program:
        s_common(parser)
        parser.add_argument("--fq", help="R2 read fastq.", required=True)


class Mapping_tag(Mt):
    pass


def mapping_tag(args):
    mt(args)
