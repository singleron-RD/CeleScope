from celescope.vdj.__init__ import CHAINS
from celescope.tools.report import reporter
from celescope.tools.utils import format_number, gen_stat, log
import os
import logging
import gzip
import numpy as np
import pandas as pd
import matplotlib as mpl
import pysam
import re
import json
import argparse
mpl.use('Agg')
from matplotlib import pyplot as plt


@log
def summary(input_file, alignments, type, outdir, sample, assay, debug, not_consensus):
    chains = CHAINS[type]

    '''
    # out files
    UMI_unfiltered_file = f'{outdir}/{sample}_UMI_unfiltered.tsv'
    UMI_filtered1_file = f'{outdir}/{sample}_UMI_filtered1.tsv'
    UMI_filtered2_file = f'{outdir}/{sample}_UMI_filtered2.tsv'
    '''

    UMI_count_unfiltered_file = f'{outdir}/{sample}_UMI_count_unfiltered.tsv'
    UMI_count_filtered1_file = f'{outdir}/{sample}_UMI_count_filtered1.tsv'

    stat_prefix = 'UMIs'
    if not_consensus:
        stat_prefix = 'Reads'

    # read input_file
    with pysam.FastxFile(input_file) as fh:
        index = 0
        read_row_list = []
        for entry in fh:
            attr = entry.name.split("_")
            barcode = attr[0]
            umi = attr[1]
            dic = {"readId": index, "barcode": barcode, "UMI": umi}
            read_row_list.append(dic)
            index += 1
        df_read = pd.DataFrame(read_row_list, columns=["readId", "barcode", "UMI"])
        summary.logger.info(f"{input_file} to dataframe done.")
        total_read = df_read.shape[0]

    # init row list
    mapping_summary_row_list = []

    # mapped
    alignment = pd.read_csv(alignments, sep="\t")
    alignment.readId = alignment.readId.astype(int)
    align_read = alignment.shape[0]
    df_read.readId = df_read.readId.astype(int)
    df_align = pd.merge(df_read, alignment, on="readId", how="right")

    mapping_summary_row_list.append({
        "item": f"{stat_prefix} Mapped to Any VDJ Gene",
        "count": align_read,
        "total_count": total_read,
    })

    # CDR3
    df_CDR3 = df_align[~pd.isnull(df_align["aaSeqCDR3"])]
    align_read_with_CDR3 = df_CDR3.shape[0]
    mapping_summary_row_list.append({
        "item": f"{stat_prefix} with CDR3",
        "count": align_read_with_CDR3,
        "total_count": total_read,
    })

    # correct CDR3
    df_correct_CDR3 = df_CDR3[~(df_CDR3["aaSeqCDR3"].str.contains(r"\*"))]
    align_read_with_correct_CDR3 = df_correct_CDR3.shape[0]
    mapping_summary_row_list.append({
        "item": f"{stat_prefix} with Correct CDR3",
        "count": align_read_with_correct_CDR3,
        "total_count": total_read,
    })

    # VDJ
    df_VJ = df_correct_CDR3[
        (~pd.isnull(df_correct_CDR3['bestVGene'])) &
        (~pd.isnull(df_correct_CDR3['bestJGene']))
    ]
    df_VJ = df_VJ[df_VJ.bestVGene.str[:3] == df_VJ.bestJGene.str[:3]]
    df_VJ["chain"] = df_VJ.bestVGene.str[:3]
    df_VJ["VJ_pair"] = df_VJ["bestVGene"] + "_" + df_VJ["bestJGene"]
    Reads_Mapped_Confidently_to_VJ_Gene = df_VJ.shape[0]
    mapping_summary_row_list.append({
        "item": f"{stat_prefix} Mapped Confidently to VJ Gene",
        "count": Reads_Mapped_Confidently_to_VJ_Gene,
        "total_count": total_read
    })

    # chain
    for chain in chains:
        df_chain = df_VJ[df_VJ.chain == chain]
        Reads_Mapped_to_chain = df_chain.shape[0]
        mapping_summary_row_list.append({
            "item": f"{stat_prefix} Mapped to {chain}",
            "count": Reads_Mapped_to_chain,
            "total_count": total_read,
        })

    # unique UMI
    df_UMI = df_VJ.drop_duplicates(subset=["barcode", "UMI"], keep="first")

    # filter1: keep top 1 in each combinations
    groupby_elements = [
        'barcode',
        'chain',
        'bestVGene',
        'bestJGene',
        'aaSeqCDR3',
        'nSeqCDR3',
    ]
    df_UMI_count = df_UMI.groupby(
        groupby_elements, as_index=False).agg({"UMI": "count"})
    df_UMI_count = df_UMI_count.sort_values("UMI", ascending=False)
    # out unfiltered
    df_UMI_count.to_csv(UMI_count_unfiltered_file, sep="\t", index=False)

    df_UMI_count_filter1 = df_UMI_count.groupby(
        ["barcode", "chain"], as_index=False).head(1)
    # out filtered1
    df_UMI_count_filter1.to_csv(
        UMI_count_filtered1_file,
        sep="\t",
        index=False)

    if debug:
        unique_UMI = df_UMI.shape[0]
        mapping_summary_row_list.append({
            "item": "UMI unique count",
            "count": unique_UMI,
            "total_count": align_read_with_correct_CDR3,
        })
        UMI_after_Contamination_Filtering = df_UMI_count.filter1.UMI.sum()
        mapping_summary_row_list.append({
            "item": "UMI after Contamination Filtering",
            "count": UMI_after_Contamination_Filtering,
            "total_count": unique_UMI,
        })

    # stat file
    df = pd.DataFrame(
        mapping_summary_row_list,
        columns=[
            "item",
            "count",
            "total_count"])
    stat_file = f'{outdir}/stat.txt'
    gen_stat(df, stat_file)

    # report
    STEP = 'mapping_vdj'
    name = f'{type}_{STEP}'
    t = reporter(
        name=name,
        sample=sample,
        stat_file=stat_file,
        outdir=outdir + '/..',
        assay=assay,
    )
    t.get_report()


@log 
def mixcr(outdir, sample, input_file, thread, species):
    report = f"{outdir}/{sample}_align.txt"
    not_align_fq = f"{outdir}/not_align.fq"
    read2_vdjca = f"{outdir}/read2.vdjca"
    alignments = f"{outdir}/{sample}_alignments.txt"

    cmd = f"""
mixcr align \
--force-overwrite \
--species {species} \
-t {thread} \
--not-aligned-R1 {not_align_fq} \
--report {report} \
-OallowPartialAlignments=true \
-OvParameters.geneFeatureToAlign=VTranscriptWithP \
{input_file} \
{read2_vdjca}
mixcr exportAlignments \
{read2_vdjca} {alignments} \
-readIds --force-overwrite -vGene -dGene -jGene -cGene \
-nFeature CDR3 -aaFeature CDR3\n"""
    mixcr.logger.info(cmd)
    os.system(cmd)
    return alignments


@log
def mapping_vdj(args):
    sample = args.sample
    outdir = args.outdir
    fq = args.fq
    type = args.type
    debug = args.debug
    assay = args.assay
    thread = int(args.thread)
    not_consensus = args.not_consensus
    species = args.species

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % outdir)

    input_file = fq
    alignments = mixcr(outdir, sample, input_file, thread, species)

    # summary
    summary(input_file, alignments, type, outdir, sample, assay, debug, not_consensus)


def get_opts_mapping_vdj(parser, sub_program):
    parser.add_argument("--type", help='TCR or BCR', required=True)
    parser.add_argument("--debug", action='store_true')
    parser.add_argument('--species', choices=['hs', 'mmu'], help='human or mouse', default='hs')
    parser.add_argument("--not_consensus", action='store_true', help="input fastq is not consensus")
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--fq", required=True)
        parser.add_argument("--thread", default=1)
        parser.add_argument('--assay', help='assay', required=True)
