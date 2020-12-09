from celescope.fusion.__init__ import __STEPS__, __ASSAY__


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

    step = "STAR_fusion"
    args.input_read = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.fusion.STAR_fusion import STAR_fusion
    STAR_fusion(args)

    step = 'count_fusion'
    args.outdir = f'{outdir_dic[step]}/'
    args.bam = f'{outdir_dic["STAR_fusion"]}/{sample}_Aligned.sortedByCoord.out.bam'
    from celescope.fusion.count_fusion import count_fusion
    count_fusion(args)
