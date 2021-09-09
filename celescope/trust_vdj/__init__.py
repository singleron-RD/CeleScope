__STEPS__ = ['sample', 'barcode', 'cutadapt', 'assemble', 'summarize']
__ASSAY__ = 'trust_vdj'

TOOLS_DIR = '/SGRNJ03/randd/zhouxin/software/TRUST4/scripts'
INDEX = '/SGRNJ03/randd/zhouxin/data_base/index'
CHAIN = {
	'TCR': ['TRA', 'TRB'], 
	'BCR': ['IGH', 'IGL', 'IGK']
	}
PAIRED_CHAIN = {
	'TCR': ['TRA_TRB'],
	'BCR': ['IGH_IGL', 'IGH_IGK']
}
CONDA_PATH = '/SGRNJ/Public/Software/conda_env/zhouxinT'