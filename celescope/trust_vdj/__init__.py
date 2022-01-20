import os
import celescope.tools

__STEPS__ = ['sample', 'barcode', 'cutadapt', 'assemble', 'summarize', 'mapping']
__ASSAY__ = 'trust_vdj'

TOOLS_DIR = os.path.dirname(celescope.tools.__file__) + '/trust4'
INDEX = f'{TOOLS_DIR}/database'

TOOLS_DIR = '/SGRNJ03/randd/cjj/zhouxin/zhouxin/software/TRUST4/scripts'
INDEX = '/SGRNJ03/randd/cjj/zhouxin/zhouxin/data_base/index'
CHAIN = {
	'TCR': ['TRA', 'TRB'], 
	'BCR': ['IGH', 'IGL', 'IGK']
	}
PAIRED_CHAIN = {
	'TCR': ['TRA_TRB'],
	'BCR': ['IGH_IGL', 'IGH_IGK']
}
CONDA_PATH = '/SGRNJ/Public/Software/conda_env/vdj_trust4'