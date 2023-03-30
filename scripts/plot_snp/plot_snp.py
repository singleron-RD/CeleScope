from celescope.snp.analysis_snp import parse_variant_ann
from celescope.tools.target_metrics import get_gene_list


def keep_in_gene(snpeff_ann_vcf_file, final_vcf_file, gene_list):
    """
    Output:
        final_vcf_file
    """
    gene_list, _, _ = parse_variant_ann(snpeff_ann_vcf_file)
    with pysam.VariantFile(snpeff_ann_vcf_file) as vcf_in:
        with pysam.VariantFile(final_vcf_file, 'w', header=vcf_in.header) as vcf_out:
            for i, record in enumerate(vcf_in.fetch()):
                if gene_list[i] in gene_list:
                    vcf_out.write(record)

