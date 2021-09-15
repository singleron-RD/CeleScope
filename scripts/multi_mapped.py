import itertools
import collections
import sys

import pysam

name_sorted_bam = sys.argv[1]

gene_counts = collections.defaultdict(int)
with pysam.AlignmentFile(name_sorted_bam, "rb") as bam:
    for _name, group in itertools.groupby(bam, lambda read: read.query_name):
        if len(group) == 1:
            continue
        genes = set()
        for read in group:
            if not read.has_tag("GN"):
                continue
            for gene in read.get_tag("GN").split(","):
                genes.add(gene)
        gene_num = len(genes)
        for gene in genes:
            gene_counts[gene] += 1 / gene_num

sorted_counts = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
for item in sorted_counts:
    print(item)
