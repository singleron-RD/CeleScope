import re

from celescope.tools import utils
from celescope.tools.matrix import Features

class GtfParser:
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file
        self.gene_id = []
        self.gene_name = []
        self.id_name = {}

        self.load_gtf()

    @utils.add_log
    def load_gtf(self):
        """
        get gene_id:gene_name from gtf file
            - one gene_name with multiple gene_id: "_{count}" will be added to gene_name.
            - one gene_id with multiple gene_name: error.
            - duplicated (gene_name, gene_id): ignore duplicated records and print a warning.
            - no gene_name: gene_id will be used as gene_name.

        Returns:
            {gene_id: gene_name} dict
        """

        gene_id_pattern = re.compile(r'gene_id "(\S+)";')
        gene_name_pattern = re.compile(r'gene_name "(\S+)"')

        with utils.generic_open(self.gtf_file) as f:
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('#'):
                    continue
                tabs = line.split('\t')
                gtf_type, attributes = tabs[2], tabs[-1]
                # TODO: fix exon without gene_id
                if gtf_type == 'gene':
                    gene_id = gene_id_pattern.findall(attributes)[-1]
                    gene_names = gene_name_pattern.findall(attributes)
                    if not gene_names:
                        gene_name = gene_id
                    else:
                        gene_name = gene_names[-1]

                    if gene_id in self.id_name:
                        assert self.id_name[gene_id] == gene_name, (
                            'one gene_id with multiple gene_name '
                            f'gene_id: {gene_id}, '
                            f'gene_name this line: {gene_name}'
                            f'gene_name previous line: {self.id_name[gene_id]}'
                        )
                        self.load_gtf.logger.warning(
                            'duplicated (gene_id, gene_name)'
                            f'gene_id: {gene_id}, '
                            f'gene_name {gene_name}'
                        )
                    else:
                        self.gene_id.append(gene_id)
                        self.gene_name.append(gene_name)
                        self.id_name[gene_id] = gene_name
                    
    def get_features(self):
        """
        Returns:
            Features object
        """
        features = Features(self.gene_id, self.gene_name)
        return features

