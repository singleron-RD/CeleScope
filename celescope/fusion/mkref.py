from celescope.tools.make_ref import MakeRef_STAR


class Mkref_fusion(MakeRef_STAR):
    """
    ## Features
    - Create a fusion genome directory.

    ## Output

    - STAR genome index files
    - Genome config file

    ## Usage
    ```
    celescope fusion mkref \\
    --genome_name {genome_name} \\
    --fasta fusion.fasta \\
    --fusion_pos fusion_pos.txt \\
    ```
    """

    def __init__(self, genome_type, args):
        super().__init__(genome_type, args)
        self.files["fusion_pos"] = args.fusion_pos


def mkref(args):
    genome_type = "fusion"
    with Mkref_fusion(genome_type, args) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    MakeRef_STAR.opts(parser, sub_program)
    if sub_program:
        parser.add_argument(
            "--fusion_pos",
            help="""
fusion position file. A two column tab-delimited text file with header.
"pos" is the end postion of the first gene(1-based).
e.g.  
name\tpos  
PML_3\t183  
PML_4\t254  
PML_5\t326  
PML_6\t204   
""",
            required=True,
        )
