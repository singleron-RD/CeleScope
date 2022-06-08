import subprocess
import pandas as pd

from celescope.tools import utils
from celescope.flv_trust4_split.__init__ import INDEX, TOOLS_DIR


"""
trust goes through the following steps:
			0: start from beginning (candidate read extraction)
			1: start from assembly
			2: start from annotation
			3: start from generating the report table
"""
@utils.add_log
def extract_candidate_reads(species, index_prefix, outdir, sample, fq1, fq2, barcodeRange, umiRange, single_thread):
    """
    extract reads map to index_prefix
    index_prefix can be 'bcrtcr' + chains
    """
    cmd = (
        f'fastq-extractor -t {single_thread} '
        f'-f {INDEX}/{species}/{index_prefix}.fa '
        f'-o {outdir}/{sample}_{index_prefix} '
        f'--barcodeStart {barcodeRange[0]} '
        f'--barcodeEnd {barcodeRange[1]} '
        f'--umiStart {umiRange[0]} '
        f'--umiEnd {umiRange[1]} '
        f'-u {fq2} '
        f'--barcode {fq1} '
        f'--UMI {fq1} '
        )
    extract_candidate_reads.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def get_trust_report(sample, filedir):
    cmd = (
        f'perl {TOOLS_DIR}/trust-simplerep.pl '
        f'{filedir}/{sample}_cdr3.out > '
        f'{filedir}/{sample}_report.tsv'
    )
    get_trust_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def filter_trust_report(sample ,filedir):
    cmd = f''' awk '$4!~"_" && $4!~"?"' {filedir}/{sample}_report.tsv > {filedir}/{sample}_filter_report.tsv '''
    filter_trust_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def get_bc_report(sample, filedir):
    cmd = (
        f'perl {TOOLS_DIR}/trust-barcoderep.pl '
        f'{filedir}/{sample}_cdr3.out > '
        f'{filedir}/{sample}_barcode_report.tsv ' 
    )
    get_bc_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def get_bcfilter_report(sample, filedir):
    cmd = (
        f'python {TOOLS_DIR}/barcoderep-filter.py '
        f'-b {filedir}/{sample}_barcode_report.tsv > '
        f'{filedir}/{sample}_barcode_filter_report.tsv '
    )
    get_bcfilter_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def convert_barcode_report(barcode_report, outdir):
    cmd = (
        f'{TOOLS_DIR}/trust-barcoderep-to-10X.pl '
        f'{barcode_report} '
        f'{outdir} '
    )
    convert_barcode_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def trust_assemble(species, outdir, sample, single_thread, trimLevel=1):

    cmd = (
        f'trust4 -t {single_thread} '
        f'-f {INDEX}/{species}/bcrtcr.fa '
        f'-o {outdir}/{sample} '
        f'-u {outdir}/{sample}.fq '
        f'--barcode {outdir}/{sample}_bc.fa '
        f'--UMI {outdir}/{sample}_umi.fa '
        f'--trimLevel {trimLevel}'
    )
    trust_assemble.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def annotate(sample, outdir, species, single_thread):

    cmd = (
        f'annotator -f {INDEX}/{species}/IMGT+C.fa '
        f'-a {outdir}/{sample}_final.out '
        f'-t {single_thread} '
        f'-o {outdir}/{sample} '
        f'--barcode --UMI --noImpute '
        f'--readAssignment {outdir}/{sample}_assign.out '
        f'-r {outdir}/{sample}_assembled_reads.fa > {outdir}/{sample}_annot.fa'
    )
    annotate.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


def get_vj_annot(df, chains, pairs):
    fl_pro_pair_df = pd.DataFrame(df[df['productive']==True].barcode.value_counts())
    fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']>=2]
    Result = []
    cell_nums = len(set(df['barcode'].tolist()))
    Result.append({
        'name': 'Cells With Productive V-J Spanning Pair',
        'value': fl_pro_pair_df.shape[0],
        'total': cell_nums,
    })
    for p in pairs:
        chain1 = p.split('_')[0]
        chain2 = p.split('_')[1]
        cbs1 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain1)].barcode.tolist())
        cbs2 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain2)].barcode.tolist())
        paired_cbs = len(cbs1.intersection(cbs2))
        Result.append({
            'name': f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
            'value': paired_cbs,
            'total': cell_nums,
            'help_info': "Fraction of cell-associated barcodes with one productive contig for each chain of the receptor pair.A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
        })
    for c in chains:
        Result.append({
            'name': f'Cells With {c} Contig',
            'value': len(set(df[df['chain']==c].barcode.tolist())),
            'total': cell_nums,
            'help_info': f'Fraction of cell-associated barcodes with at least one {c} contig annotated as a full or partial V(D)J gene'
        })
        Result.append({
            'name': f'Cells With CDR3-annotated {c} Contig',
            'value': len(set(df[(df['chain']==c)&(df['cdr3']!=None)].barcode.tolist())),
            'total': cell_nums,
        })
        Result.append({
            'name': f'Cells With V-J Spanning {c} Contig',
            'value': len(set(df[(df['full_length']==True)&(df['chain']==c)].barcode.tolist())),
            'total': cell_nums,
            'help_info': f"Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for {c}"
        })
        Result.append({
            'name': f'Cells With Productive {c} Contig',
            'value': len(set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==c)].barcode.tolist())),
            'total': cell_nums,
            'help_info': "Fraction of cell-associated barcodes with productive IGL chain. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
        })

    return Result