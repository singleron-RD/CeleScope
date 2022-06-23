import os
import celescope.tools


__STEPS__ = ['sample', 'barcode', 'mapping', 'assemble', 'summarize', 'annotation']
__ASSAY__ = 'flv_trust4'

TOOLS_DIR = os.path.dirname(celescope.tools.__file__) + '/trust4'
REF_DIR = f'{TOOLS_DIR}/database'

CHAIN = {
	'TCR': ['TRA', 'TRB'], 
	'BCR': ['IGH', 'IGL', 'IGK']
	}
PAIRED_CHAIN = {
	'TCR': ['TRA_TRB'],
	'BCR': ['IGK_IGH', 'IGL_IGH']
}