from celescope.tools import utils
from celescope.tools import reference


class Mkgtf:
    """
    ##Features

    - Filter GTF files.

    GTF files downloaded from sites like ENSEMBL and UCSC often contain transcripts and genes which need to be filtered.
    Celescope provides mkgtf, a simple utility to filter genes. The command syntax requires input and output fitlered GTF file names.

    The following gene biotypes will be kept.
    ```
    protein_coding
    lncRNA
    antisense
    IG_LV_gene
    IG_V_gene
    IG_V_pseudogene
    IG_D_gene
    IG_J_gene
    IG_J_pseudogene
    IG_C_gene
    IG_C_pseudogene
    TR_V_gene
    TR_V_pseudogene
    TR_D_gene
    TR_J_gene
    TR_J_pseudogene
    TR_C_gene
    ```
    """

    def __init__(self, in_gtf, out_gtf, attributes, add_intron=True):
        self.in_gtf_fn = in_gtf
        self.out_gtf_fn = out_gtf
        self.add_intron = add_intron

        self.attributes = {}
        for attr_str in attributes.split(";"):
            if attr_str:
                attr, val = attr_str.split("=")
                val = set(val.split(","))
                self.attributes[attr] = val

    @utils.add_log
    def run(self):
        runner = reference.GtfBuilder(
            self.in_gtf_fn, self.out_gtf_fn, self.attributes, add_intron=self.add_intron
        )
        runner.build_gtf()


def mkgtf(args):
    runner = Mkgtf(args.gtf, args.out_gtf, args.attributes, not args.skip_intron)
    runner.run()


def get_opts_mkgtf(parser, sub_program=True):
    if sub_program:
        parser.add_argument("gtf", help="raw gtf file")
        parser.add_argument("out_gtf", help="output gtf file")
        parser.add_argument(
            "--attributes",
            help="Attributes to keep. Example: `gene_biotype=protein_coding,lncRNA,antisense;`",
            default="gene_biotype=protein_coding,lncRNA,antisense,IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene;",
        )
        parser.add_argument(
            "--skip_intron", action="store_true", help="Do not add intron to gtf"
        )

    return parser
