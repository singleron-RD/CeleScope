from celescope.vdj.__init__ import __STEPS__, __ASSAY__


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

    step = 'mapping_vdj'
    args.fq = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.vdj.mapping_vdj import mapping_vdj
    mapping_vdj(args)

    step = 'count_vdj'
    args.UMI_count_filter1_file = f'{outdir_dic["mapping_vdj"]}/{sample}_UMI_count_filtered1.tsv'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.vdj.count_vdj import count_vdj
    count_vdj(args)
