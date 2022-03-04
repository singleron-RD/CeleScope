import os
from collections import OrderedDict

__VERSION__ = "1.8.0b0"
__version__ = __VERSION__

ASSAY_LIST = [
    "rna",
    'vdj',
    'tag',
    'dynaseq',
    'snp',
    'capture_virus',
    'fusion',
    'hla',
    'capture_rna',
    'citeseq',
    'tcr_fl',
]

ROOT_PATH = os.path.dirname(__file__)

RELEASED_ASSAYS = ['rna', 'vdj', 'tag', 'dynaseq', 'snp', 'capture_virus']

# argument help
HELP_DICT = {
    'match_dir': 'Match celescope scRNA-Seq directory.',
    'gene_list': 'Required. Gene list file, one gene symbol per line. Only results of these genes are reported.',
    'genomeDir': 'Required. Genome directory after running `celescope rna mkref`.',
    'thread': 'Thread to use.',
    'debug': 'If this argument is used, celescope may output addtional file for debugging.',
    'fasta': 'Required. Genome fasta file. Use absolute path or relative path to `genomeDir`.',
    'outdir': 'Output directory.',
    'matrix_dir': 'Match celescope scRNA-Seq matrix directory.',
    'panel': 'The prefix of bed file in `celescope/data/snp/panel/`, such as `lung_1`.',
    'virus_genomeDir': 'Required. Virus genome directory after running `celescope capture_virus mkref`.',
    'threshold_method': 'One of [otsu, auto, hard, none].',
    'tsne_file': 'match_dir t-SNE coord file. Do not required when `--match_dir` is provided.',
    'df_marker_file': 'match_dir df_marker_file. Not required when `--match_dir` is provided.',
    'citeseq_mtx':'citeseq matrix .gz file',
    'strand':'gene strand file, the format is "geneID,+/-"',
    'bam_for_conversion':'featureCount bam(sortedByCoord), must have "MD" tag, set in star step',
    'cell':'barcode cell list',
    'tsne_for_replace_tsne':'tsne file from analysis step',
    'mat_for_replace_tsne':'matrix replacement file, from replacement step',
    'rep_for_replace_tsne':'cell replacement file, from replacement step',
    'mincell_for_replace_tsne':'turn-over in at least cells, default 5',
    'topgene_for_replace_tsne':'show top N genes,default 10',
    'bg_cov_for_replacement':'background snp depth filter, lower than bg_cov will be discarded. Only valid in csv format',
    'bam_for_replacement':'bam file from conversion step',
    'bg_for_replacement':'background snp file, csv or vcf format',
    'cell_keep_for_replacement':'filter cell',
    'min_cell_for_replacement':'a gene expressed in at least cells, default 10',
    'min_gene_for_replacement':'at least gene num in a cell, default 10',
    'fusion_genomeDir':'Fusion genome directory.',
    'flanking_base':'Number of bases flanking the fusion position.',
    'fusion_pos':"""
fusion position file. A two column tab-delimited text file with header.
"pos" is the end postion of the first gene(1-based).
e.g.  
name\tpos  
PML_3\t183  
PML_4\t254  
PML_5\t326  
PML_6\t204   
""",
    'genomeSAindexNbases':'STAR genomeSAindexNbases',
    'gtf':'Required. Genome gtf file. Use absolute path or relative path to `genomeDir`.',
    'mt_gene_list':"""Mitochondria gene list file. Use absolute path or relative path to `genomeDir`.
It is a plain text file with one gene per line. 
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.""",
    'virus_file':'virus UMI count file',
    'hard_threshold':'int, use together with `--threshold_method hard`',
    'vcf':'norm vcf file',
    'annovar_config':'ANNOVAR config file.',
    'split_fastq':'If used, will split scRNA-Seq fastq file according to tag assignment.',
    'split_matrix':'If used, will split scRNA-Seq matrix file according to tag assignment.',
    'split_vdj':'If used, will split scRNA-Seq vdj count file according to tag assignment.',
    'vdj_dir':'Match celescope vdj directory. Required when --split_vdj is specified.',
    'umi_tag_file':'UMI tag file.',
    'R1_read':'R1 read path.',
    'UMI_count_filter_file':'Required. File from step mapping_vdj.',
    'type':'Required. `TCR` or `BCR`. ',
    'UMI_min':'Default `auto`. Minimum UMI number to filter. The barcode with UMI>=UMI_min is considered to be cell.',
    'species':'Default `hs`. `hs`(human) or `mmu`(mouse). ',
    'fq':'Required. Input fastq file.',
    'not_consensus':'Input fastq is not consensused.',
}

# report metrics help
HELP_INFO_DICT = {
    'matched_barcode_number': {
        'display': 'Number of Matched Cells',
        'info': 'cell barcode number of matched scRNA-Seq sample',
    }
}
