from celescope.tools.tag.mapping_tag import Mapping_tag
from celescope.tools.step import s_common

def get_opts_mapping_tag(parser, sub_program):
    parser.add_argument(
        "--fq_pattern",
        help="""R2 read pattern. The number after the letter represents the number of bases. The `fq_pattern` of CLindex is `L25C15`
`L` linker(common sequences)  
`C` tag barcode  
""",
        default='L25C15'
    )
    parser.add_argument(
        "--barcode_fasta",
        help="""Required. Tag barcode fasta file. It will check the mismatches between tag barcode 
sequence in R2 reads with all tag barcode sequence in barcode_fasta. 
It will assign read to the tag with mismatch < threshold. 
If no such tag exists, the read is classified as invalid.

You can find the barcode fasta file under `celescope/data/Clindex`
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

def mapping_tag(args):
    with Mapping_tag_tag(args) as runner:
        runner.run()


class Mapping_tag_tag(Mapping_tag):
    '''
    ## Features
    - Assign tag to R2 reads.
    ## Output

    - `{sample}_invalid_barcode.tsv` Reads count with invalid tag.

    - `{sample}_read_count.tsv` Reads count with effective tag in each barcode.
    '''