import re
import collections
import csv
import sys

from celescope.tools import utils
from celescope.tools.matrix import Features


class GeneIdNotFound(Exception):
    pass

class GtfParser:
    def __init__(self, gtf_fn):
        self.gtf_fn = gtf_fn
        self.gene_id = []
        self.gene_name = []
        self.id_name = {}

    def get_properties_dict(self, properties_str):
        
        if isinstance(properties_str, dict):
            return properties_str

        properties = collections.OrderedDict()
        pattern = re.compile(r'(\S+?)\s*"(.*?)"')
        for m in re.finditer(pattern, properties_str):
            key = m.group(1)
            value = m.group(2)
            properties[key] = value

        return properties

    def gtf_reader_iter(self):
        """
        Yield:
            row: str
            is_comment: bool
            annotation: str; exon, gene, etc
            properties: dict; gtf properties in the last column(9th col)
        """
        with utils.generic_open(self.gtf_fn, mode='rt') as f:
            reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(reader, start=1):
                if len(row) == 0:
                    continue
                if row[0].startswith('#'):
                    yield row, True, None, None
                    continue

                if len(row) != 9:
                    sys.exit(f"Invalid number of columns in GTF line {i}: {row}\n")

                if row[6] not in ['+', '-']:
                    sys.exit(f"Invalid strand in GTF line {i}: {row}\n")
                
                properties = self.get_properties_dict(row[8])
                annotation = row[2]
                if annotation == 'exon':
                    if 'gene_id' not in properties:
                        raise GeneIdNotFound(f"Property 'gene_id' not found in GTF line {i}: {row}\n")

                yield row, False, annotation, properties

    @utils.add_log
    def get_id_name(self):
        """
        get gene_id:gene_name from gtf file
            - one gene_name with multiple gene_id: allowed.
            - one gene_id with multiple gene_name: error.
            - duplicated (gene_name, gene_id): ignore duplicated records and print a warning.
            - no gene_name: gene_id will be used as gene_name.

        """
        for _row, _is_comment, annotation, properties in self.gtf_reader_iter():
            if annotation == 'gene':
                gene_id = properties['gene_id']
                if 'gene_name' not in properties:
                    gene_name = gene_id
                else:
                    gene_name = properties['gene_name']

                if gene_id in self.id_name:
                    assert self.id_name[gene_id] == gene_name, (
                        'one gene_id with multiple gene_name '
                        f'gene_id: {gene_id}, '
                        f'gene_name this line: {gene_name}'
                        f'gene_name previous line: {self.id_name[gene_id]}'
                    )
                    self.get_id_name.logger.warning(
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
        if not self.gene_id:
            sys.exit("Empty self.gene_id. Run self.get_id_name first.")
        features = Features(self.gene_id, self.gene_name)
        return features


class GtfBuilder:
    def __init__(self, in_gtf_fn, out_gtf_fn, attributes):
        self.in_gtf_fn = in_gtf_fn
        self.out_gtf_fn = out_gtf_fn
        self.attributes = attributes

    @utils.add_log
    def build_gtf(self):
        self.build_gtf.logger.info("Writing filtered GTF file...")
        gp = GtfParser(self.in_gtf_fn)
        with open(self.out_gtf_fn, 'w') as f:
            writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
            for row, is_comment, _, properties in gp.gtf_reader_iter():
                if is_comment:
                    writer.writerow(row)
                    continue

                remove = False
                for key, value in properties.items():
                    if key in self.attributes and value not in self.attributes[key]:
                        remove = True

                if not remove:
                    writer.writerow(row)