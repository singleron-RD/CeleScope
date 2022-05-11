import os
import celescope.tools


__STEPS__ = ['sample', 'barcode', 'cutadapt', 'assemble', 'summarize', 'annotation']
__ASSAY__ = 'fl_vdj_TRUST4'

TOOLS_DIR = os.path.dirname(celescope.tools.__file__) + '/trust4'
INDEX = f'{TOOLS_DIR}/database'

CHAIN = {
	'TCR': ['TRA', 'TRB'], 
	'BCR': ['IGH', 'IGL', 'IGK']
	}
PAIRED_CHAIN = {
	'TCR': ['TRA_TRB'],
	'BCR': ['IGH_IGL', 'IGH_IGK']
}