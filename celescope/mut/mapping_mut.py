import os
from celescope.tools.utils import format_number, log


@log
def mapping_mut(args):

    sample = args.sample
    outdir = args.outdir
    fq = args.fq
    indel_genomeDir = args.indel_genomeDir
    runThreadN = args.thread

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    out_prefix = f'{outdir}/{sample}_'
    out_BAM = out_prefix + "Aligned.sortedByCoord.out.bam"

    cmd = f"STAR \
 --genomeDir {indel_genomeDir} \
 --readFilesIn {fq}\
 --readFilesCommand zcat\
 --outSAMtype BAM SortedByCoordinate\
 --runThreadN {runThreadN}\
 --limitBAMsortRAM 10000000000\
 --outFileNamePrefix {out_prefix}"

    mapping_mut.logger.info(cmd)
    os.system(cmd)

    # samtools index
    cmd = "samtools index {out_BAM}".format(out_BAM=out_BAM)
    mapping_mut.logger.info(cmd)
    os.system(cmd)


def get_opts_mapping_mut(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--fq", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--indel_genomeDir',
        help='insertion or deletion STAR indexed genome directory',
        required=True)
    parser.add_argument("--thread", help='STAR thread', default=1)