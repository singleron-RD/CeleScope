import os
import glob

__STEPS__ = ['sample', 'convert', 'vdj_10X', ]
__ASSAY__ = 'vdj10X'

soft_version = ['3.0.2', '3.1.0', '4.0.0', '6.0.0']
default_version = '4.0.0'
species_choice = ['hs', 'mmu']
soft_dict = {}
ref_dict = {}
for version in soft_version:
    soft_dict[version] = f'/SGRNJ/Database/script/soft/cellranger/cellranger-{version}/cellranger'
    ref_dict[version] = {}
    for species in species_choice:
        default_ref = glob.glob(
            f'/SGRNJ/Database/script/soft/cellranger/vdj_ref/{default_version}/{species}/refdata*{default_version}')[0]
        ref_path = glob.glob(f'/SGRNJ/Database/script/soft/cellranger/vdj_ref/{version}/{species}/refdata*{version}')
        if ref_path:
            ref_dict[version][species] = ref_path[0]
        else:
            ref_dict[version][species] = default_ref



