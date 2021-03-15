import subprocess
from celescope.tools.utils import log
import os
from celescope.vdj10X.__init__ import ref_dict

@log
def vdj_10X(args):
    sample = args.sample
    outdir = args.outdir
    thread = args.thread
    mem = args.mem
    species = args.species
    ref = ref_dict[species]

    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % args.outdir)

    cmd = (
        f'/SGRNJ/Database/script/soft/cellranger/cellranger-3.0.2/cellranger vdj '
        f'--id={sample} '
        f'--reference={ref} '
        f'--fastqs=../01.convert '
        f'--sample={sample} '
        f'--localcores={thread} '
        f'--localmem={mem} '
    )
    vdj_10X.logger.info(cmd)
    with open(f'{outdir}/{sample}_vdj_10X.sh', 'w') as f:
        f.write(cmd)
    os.chdir(outdir)
    subprocess.run(cmd, shell=True)


def get_opts_vdj_10X(parser, sub_program):
    parser.add_argument('--species', help='species', choices=['hs','mmu'], required=True)
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--thread', help='number of threads', default=4)
        parser.add_argument('--mem', help='memory (G)', default=10)
    return parser