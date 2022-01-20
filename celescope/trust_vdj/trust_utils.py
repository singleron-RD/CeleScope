import subprocess
import pandas as pd
import pysam

from celescope.tools import utils
from celescope.trust_vdj.__init__ import CHAIN, INDEX, TOOLS_DIR


@utils.add_log
def VDJ_mapping(thread, species, index_prefix, outdir, sample, fq1, fq2, barcodeRange, umiRange):
    cmd = (
        f'fastq-extractor -t {thread} '
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
    VDJ_mapping.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)

@utils.add_log
def get_trust_report(filedir, sample):
    cmd = (
        f'perl {TOOLS_DIR}/trust-simplerep.pl '
        f'{filedir}/{sample}_cdr3.out > '
        f'{filedir}/report.out'
    )
    get_trust_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)

@utils.add_log
def get_bc_report(filedir, sample):
    cmd = (
        f'perl {TOOLS_DIR}/trust-barcoderep.pl '
        f'{filedir}/{sample}_cdr3.out > '
        f'{filedir}/barcoderep.tsv ' 
    )
    get_bc_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)

@utils.add_log
def get_bcfilter_report(filedir, sample):
    cmd = (
        f'python {TOOLS_DIR}/barcoderep-filter.py '
        f'-b {filedir}/barcoderep.tsv > '
        f'{filedir}/barcoderepfl.tsv '
    )
    get_bcfilter_report.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def trust_assemble(thread, species, outdir, sample, trimLevel=1):

    cmd = (
        f'trust4 -t {thread} '
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
def get_full_len_assembly(filedir, sample):
    cmd = (
        f'perl {TOOLS_DIR}/GetFullLengthAssembly.pl '
        f'{filedir}/{sample}_annot.fa > '
        f'{filedir}/{sample}_full_len.fa '
    )
    get_full_len_assembly.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def annotate(sample, thread, outdir, species):

    cmd = (
        f'annotator -f {INDEX}/{species}/IMGT+C.fa '
        f'-a {outdir}/{sample}_final.out '
        f'-t {thread} '
        f'-o {outdir}/{sample} '
        f'--barcode --UMI --noImpute '
        f'--readAssignment {outdir}/{sample}_assign.out '
        f'-r {outdir}/{sample}_assembled_reads.fa > {outdir}/{sample}_annot.fa'
    )
    annotate.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def fa_to_csv(outdir, sample):
    # file name
    full_len_fa = f'{outdir}/{sample}_full_len.fa'
    assign_file = f'{outdir}/{sample}_assign.out'
    # reads assignment 
    assignment = pd.read_csv(assign_file, sep='\t', header=None)
    assignment['read_barcode'] = assignment[0].apply(lambda x: x.split('_')[0])
    assignment['contig_barcode'] = assignment[1].apply(lambda x: x.split('_')[0])
    assignment['match_barcode'] = assignment[['read_barcode', 'contig_barcode']].apply(lambda x: x['read_barcode']==x['contig_barcode'], axis=1)
    assignment = assignment[assignment['match_barcode']==True]
    assignment['umi'] = assignment[0].apply(lambda x: x.split('_')[1])
    # write contig csv
    contigs = open(f'{outdir}/{sample}_contig.csv', 'w')
    # contigs.write('barcode\tis_cell\tcontig_id\thigh_confidence\tlength\tchain\tv_gene\td_gene\tj_gene\tc_gene\tfull_length\tproductive\tcdr3\tcdr3_nt\treads\tumis\traw_clonotype_id\traw_consensus_id\n')
    process_read = 0
    with pysam.FastxFile(full_len_fa) as fa:
        for read in fa:
            name = read.name
            comment = read.comment
            attrs = comment.split(' ')
            barcode = name.split('_')[0]
            is_cell = 'True'
            high_confidence = 'True'
            length = attrs[0]
            chain = attrs[2][:3]
            full_length = 'True'
            v_gene = attrs[2].split('(')[0]
            d_gene = attrs[3]
            j_gene = attrs[4].split('(')[0]
            c_gene = attrs[5]
            cdr3 = attrs[8].split('=')[1]
            cdr3_aa = 'None'
            productive = 'False'
            temp = assignment[assignment[1]==name]
            reads = str(len(temp[0].tolist()))
            umis = str(len(set(temp['umi'].tolist())))
            #raw_consensus_id = 'None'
            #raw_clonotype_id = 'None'

            string = '\t'.join([barcode, is_cell, name, high_confidence, length, chain, v_gene, d_gene, j_gene, c_gene, full_length, productive, cdr3_aa, cdr3, reads, umis])
            contigs.write(f'{string}\n')
            process_read+=1
            if process_read % 10000 == 0:
                fa_to_csv.logger.info(f'Processed {process_read} contigs')

    contigs.close()


def get_vj_annot(df, chains, pairs):
    fl_pro_pair_df = pd.DataFrame(df[df['productive']==True].barcode.value_counts())
    fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']>=2]
    Result = []
    cell_nums = len(set(df['barcode'].tolist()))
    Result.append({
        'item': 'Cells With Productive V-J Spanning Pair',
        'count': fl_pro_pair_df.shape[0],
        'total_count': cell_nums
    })
    for p in pairs:
        chain1 = p.split('_')[0]
        chain2 = p.split('_')[1]
        cbs1 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain1)].barcode.tolist())
        cbs2 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain2)].barcode.tolist())
        paired_cbs = len(cbs1.intersection(cbs2))
        Result.append({
            'item': f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
            'count': paired_cbs,
            'total_count': cell_nums
        })
    for c in chains:
        Result.append({
            'item': f'Cells With {c} Contig',
            'count': len(set(df[df['chain']==c].barcode.tolist())),
            'total_count': cell_nums
        })
        Result.append({
            'item': f'Cells With CDR3-annotated {c} Contig',
            'count': len(set(df[(df['chain']==c)&(df['productive']==True)].barcode.tolist())),
            'total_count': cell_nums
        })
        Result.append({
            'item': f'Cells With V-J Spanning {c} Contig',
            'count': len(set(df[(df['full_length']==True)&(df['chain']==c)].barcode.tolist())),
            'total_count': cell_nums
        })
        Result.append({
            'item': f'Cells With Productive {c} Contig',
            'count': len(set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==c)].barcode.tolist())),
            'total_count': cell_nums
        })

    return Result