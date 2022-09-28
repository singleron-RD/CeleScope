from celescope.tools.tag.mapping_tag import Mapping_tag, get_opts_mapping_tag

def get_opts_mapping_tag_citeseq(parser, sub_program):
    get_opts_mapping_tag(parser, sub_program)

def mapping_tag(args):
    with Mapping_tag_citeseq(args) as runner:
        runner.run()

class Mapping_tag_citeseq(Mapping_tag):
    '''
    ## Features
    - Assign tag to R2 reads.
    ## Output

    - `{sample}_invalid_barcode.tsv` Reads count with invalid tag.

    - `{sample}_read_count.tsv` Reads count with effective tag in each barcode.
    '''