from celescope.tools.tag.count_tag import Count_tag as Ct,  count_tag as ct, get_opts_count_tag as opts




class Count_tag(Ct):
    """
    ## Features

    ## Output

    - `{sample}_umi_tag.tsv` UMI count assigned to each barcode

    - `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation.
    """

def count_tag(args):
    ct(args)

def get_opts_count_tag(parser, sub_program):
    opts(parser, sub_program)