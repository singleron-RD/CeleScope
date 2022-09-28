from celescope.tools.tag.count_tag import Count_tag as Ct,  count_tag, get_opts_count_tag

def get_opts_count_tag_sweetseq(parser, sub_program):
    get_opts_count_tag(parser, sub_program)

def count_tag_sweetseq(args):
    count_tag(args)

class Count_tag(Ct):
    """
    ## Features
    - Assign tag to each cell barcode and summarize.
    ## Output

    - `{sample}_umi_tag.tsv` UMI count assigned to each barcode

    - `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation.
    """