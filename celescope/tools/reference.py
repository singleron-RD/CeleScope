import re
import collections
import csv
import sys

import pandas as pd

from celescope.tools import utils
from celescope.tools.matrix import Features

PATTERN = re.compile(r'(\S+?)\s*"(.*?)"')
gtf_row = collections.namedtuple(
    "gtf_row", "seqname source feature start end " "score strand frame attributes"
)


def row2list(row):
    """
    row: gtf_row
    """
    attr_str = "; ".join(['%s "%s"' % attr for attr in row.attributes.items()])
    return [
        row.seqname,
        row.source,
        row.feature,
        row.start,
        row.end,
        row.score,
        row.strand,
        row.frame,
        attr_str,
    ]


class GtfParser:
    def __init__(self, gtf_fn):
        self.gtf_fn = gtf_fn
        self.gene_id = []
        self.gene_name = []
        self.id_name = {}
        self.id_strand = {}

    def get_properties_dict(self, properties_str):
        """
        allow no space after semicolon
        """

        if isinstance(properties_str, dict):
            return properties_str

        properties = collections.OrderedDict()
        attrs = properties_str.split(";")
        for attr in attrs:
            if attr:
                m = re.search(PATTERN, attr)
                if m:
                    key = m.group(1).strip()
                    value = m.group(2).strip()
                    properties[key] = value

        return properties

    def gtf_reader_iter(self):
        """
        Yield:
            row: list
            gtf_row
        """
        with utils.generic_open(self.gtf_fn, mode="rt") as f:
            reader = csv.reader(f, delimiter="\t")
            for i, row in enumerate(reader, start=1):
                if len(row) == 0:
                    continue
                if row[0].startswith("#"):
                    yield row, None
                    continue

                if len(row) != 9:
                    sys.exit(f"Invalid number of columns in GTF line {i}: {row}\n")

                if row[6] not in ["+", "-"]:
                    sys.exit(f"Invalid strand in GTF line {i}: {row}\n")

                seqname = row[0]
                source = row[1]
                feature = row[2]
                # gff/gtf is 1-based, end-inclusive
                start = int(row[3])
                end = int(row[4])
                score = row[5]
                strand = row[6]
                frame = row[7]
                attributes = self.get_properties_dict(row[8])

                yield (
                    row,
                    gtf_row(
                        seqname,
                        source,
                        feature,
                        start,
                        end,
                        score,
                        strand,
                        frame,
                        attributes,
                    ),
                )

    @utils.add_log
    def get_id_name(self):
        """
        return: {gene_id:gene_name}
        """
        for _, grow in self.gtf_reader_iter():
            if not grow:
                continue
            gene_id = grow.attributes["gene_id"]
            self.id_strand[gene_id] = grow.strand
            if "gene_name" not in grow.attributes:
                gene_name = gene_id
            else:
                gene_name = grow.attributes["gene_name"]

            if gene_id not in self.id_name:
                self.gene_id.append(gene_id)
                self.gene_name.append(gene_name)
                self.id_name[gene_id] = gene_name

        return self.id_name

    def get_features(self):
        """
        Returns:
            Features object
        """
        if not self.gene_id:
            sys.exit("Empty self.gene_id. Run self.get_id_name first.")
        features = Features(self.gene_id, self.gene_name)
        return features

    def get_strand(self):
        return self.id_strand


class GtfBuilder:
    def __init__(self, in_gtf_fn, out_gtf_fn, attributes, add_intron=False):
        self.in_gtf_fn = in_gtf_fn
        self.out_gtf_fn = out_gtf_fn
        self.attributes = attributes
        self.add_intron = add_intron
        self.mt_gene_list = "mt_gene_list.txt"

    @staticmethod
    def get_introns(exons):
        sys.stderr.write("done (%d exons).\n" % len(exons))
        # add intron
        transcripts = collections.defaultdict(list)
        for grow in exons:
            if "transcript_id" in grow.attributes:
                transcripts[grow.attributes["transcript_id"]].append(grow)
        sys.stderr.write("done (%d transcripts).\n" % len(transcripts))

        introns = []
        for transcript_id, rows in transcripts.items():
            if len(set((row.seqname, row.strand) for row in rows)) != 1:
                sys.stderr.write('Malformed transcript "%s". Skipping.' % transcript_id)

            rows.sort(key=lambda row: row.start)
            starts = [row.start for row in rows]
            ends = [row.end for row in rows]

            k = 0
            for i, j in zip(ends[:-1], starts[1:]):
                if i >= j:
                    sys.stderr.write(
                        f"Skipping {transcript_id} last exon end {i} next exon start {j}\n"
                    )
                    continue

                if i + 1 <= j - 1:
                    k += 1
                    attributes = rows[0].attributes.copy()
                    if "exon_number" in attributes:
                        del attributes["exon_number"]
                    attributes["intron_number"] = k
                    intron_row = gtf_row(
                        seqname=rows[0].seqname,
                        source=rows[0].source,
                        feature="intron",
                        start=i + 1,
                        end=j - 1,
                        score=".",
                        strand=rows[0].strand,
                        frame=".",
                        attributes=attributes,
                    )

                    introns.append(intron_row)
        sys.stderr.write("done (%d introns).\n" % len(introns))
        return introns

    @utils.add_log
    def build_gtf(self):
        """
        add intron
        Filter gene biotypes
        """
        sys.stderr.write("Writing GTF file...\n")
        gp = GtfParser(self.in_gtf_fn)
        n_filter = 0
        exons = []
        mt = set()

        with open(self.out_gtf_fn, "w") as f:
            writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE, quotechar="")
            for row, grow in gp.gtf_reader_iter():
                if not grow:
                    writer.writerow(row)
                    continue

                remove = False
                for key, value in grow.attributes.items():
                    if key in self.attributes and value not in self.attributes[key]:
                        remove = True

                if not remove:
                    if grow.seqname.upper() == "MT":
                        mt.add(grow.attributes["gene_name"])
                    writer.writerow(row)
                    if grow.feature == "exon":
                        exons.append(grow)
                else:
                    n_filter += 1
            sys.stderr.write(f"filtered line number: {n_filter}\n")

            if self.add_intron:
                introns = self.get_introns(exons)
                for intron in introns:
                    writer.writerow(row2list(intron))

            if mt:
                sys.stderr.write("Writing mito genes\n")
                pd.Series(list(mt)).to_csv(self.mt_gene_list, index=False, header=False)
