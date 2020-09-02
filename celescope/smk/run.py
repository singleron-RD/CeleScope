from celescope.smk.__init__ import __STEPS__, __ASSAY__


def run(args):
    steps = __STEPS__
    args.assay = __ASSAY__
    sample = args.sample

    outdir_dic = {}
    index = 0
    for step in steps:
        outdir = f"{sample}/{index:02d}.{step}"
        outdir_dic.update({step: outdir})
        index += 1

    step = "sample"
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.tools.sample_info import sample_info
    sample_info(args)

    step = "barcode"
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.tools.barcode import barcode
    barcode(args)

    step = "cutadapt"
    args.outdir = f'{outdir_dic[step]}/'
    args.fq = f'{outdir_dic["barcode"]}/{sample}_2.fq.gz'
    from celescope.tools.cutadapt import cutadapt
    cutadapt(args)

    step = 'mapping_smk'
    args.SMK_read2 = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.smk.mapping_smk import mapping_smk
    mapping_smk(args)

    step = 'count_smk'
    args.read_file = f'{outdir_dic["mapping_smk"]}/{sample}_read_count.tsv'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.smk.count_smk import count_smk
    count_smk(args)

    step = 'analysis_smk'
    args.tsne_tag_file = f'{outdir_dic["count_smk"]}/{sample}_tsne_tag.tsv'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.smk.analysis_smk import analysis_smk
    analysis_smk(args)