from celescope.tools.tag.count_tag import Count_tag as Ct,  count_tag, get_opts_count_tag




class Count_tag(Ct):
    """
    ## Features

    ## Output

    - `{sample}_umi_tag.tsv` UMI count assigned to each barcode

    - `{sample}_tsne_tag.tsv` it is `{sample}_umi_tag.tsv` with t-SNE coordinates, gene_counts and cluster infomation.
    """
