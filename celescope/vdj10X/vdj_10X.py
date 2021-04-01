import subprocess
from celescope.tools.utils import *
import os
from celescope.vdj10X.__init__ import ref_dict, soft_dict

@add_log
def vdj_10X(args):
    sample = args.sample
    outdir = args.outdir
    thread = args.thread
    mem = args.mem
    species = args.species
    soft = args.soft

    ref_path = ref_dict[soft][species]
    soft_path = soft_dict[soft]

    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % args.outdir)

    cmd = (
        f'{soft_path} vdj '
        f'--id={sample} '
        f'--reference={ref_path} '
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
    parser.add_argument('--soft', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'], 
        default='4.0.0')
    parser.add_argument('--mem', help='memory (G)', default=10)
    if sub_program:
        s_common(parser)
    return parser