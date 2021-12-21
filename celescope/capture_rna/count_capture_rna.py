
from collections import defaultdict
from itertools import groupby

import pandas as pd
import pysam

import celescope.tools.utils as utils
from celescope.tools.count import Count, get_opts_count


class Count_capture_rna(Count):

    def bam2table(self):
        """
        read probe file
        """
        probe_gene_count_dict = utils.genDict(dim=4, valType=int)

        samfile = pysam.AlignmentFile(self.bam, "rb")
        with open(self.count_detail_file, 'wt') as fh1:
            fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'count']) + '\n')

            def keyfunc(x): return x.query_name.split('_', 1)[0]
            for _, g in groupby(samfile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                for seg in g:
                    (barcode, umi, probe) = seg.query_name.split('_')[:3]
                    if probe != 'None':
                        probe_gene_count_dict[probe]['total'][barcode][umi] += 1
                        if seg.has_tag('XT'):
                            geneID = seg.get_tag('XT')
                            geneName = self.gtf_dict[geneID]
                            probe_gene_count_dict[probe][geneName][barcode][umi] += 1
                        else:
                            probe_gene_count_dict[probe]['None'][barcode][umi] += 1
                    if not seg.has_tag('XT'):
                        continue
                    geneID = seg.get_tag('XT')
                    gene_umi_dict[geneID][umi] += 1
                for gene_id in gene_umi_dict:
                    Count.correct_umi(gene_umi_dict[gene_id])

                # output
                for gene_id in gene_umi_dict:
                    for umi in gene_umi_dict[gene_id]:
                        fh1.write('%s\t%s\t%s\t%s\n' % (barcode, gene_id, umi,
                                                        gene_umi_dict[gene_id][umi]))

        # out probe
        row_list = []
        for probe in probe_gene_count_dict:
            for geneName in probe_gene_count_dict[probe]:
                barcode_count = len(probe_gene_count_dict[probe][geneName])
                umi_count = 0
                read_count = 0
                for barcode in probe_gene_count_dict[probe][geneName]:
                    for umi in probe_gene_count_dict[probe][geneName][barcode]:
                        umi_count += len(probe_gene_count_dict[probe][geneName][barcode])
                        read_count += probe_gene_count_dict[probe][geneName][barcode][umi]
                row_list.append({
                    'probe': probe,
                    'gene': geneName,
                    'barcode_count': barcode_count,
                    'read_count': read_count,
                    'UMI_count': umi_count
                })

        df_probe = pd.DataFrame(row_list,
                                columns=['probe', 'gene', 'barcode_count', 'read_count', 'UMI_count'])
        df_probe = df_probe.groupby(['probe']).apply(
            lambda x: x.sort_values('UMI_count', ascending=False)
        )
        return df_probe

    def run(self):
        df_probe = self.bam2table()
        df_probe.to_csv(f'{self.outdir}/{self.sample}_probe_gene_count.tsv', sep='\t', index=False)

        df = pd.read_table(self.count_detail_file, header=0)
        # df_sum
        df_sum = Count.get_df_sum(df)

        # export all matrix
        self.write_matrix_10X(df, self.raw_matrix_10X_dir)

        # call cells
        cell_bc, _threshold = self.cell_calling(df_sum)

        # get cell stats
        CB_describe = self.get_cell_stats(df_sum, cell_bc)

        # export cell matrix
        df_cell = df.loc[df['Barcode'].isin(cell_bc), :]
        self.write_matrix_10X(df_cell, self.cell_matrix_10X_dir)
        (CB_total_Genes, CB_reads_count, reads_mapped_to_transcriptome) = self.cell_summary(
            df, cell_bc)

        # downsampling
        cell_bc = set(cell_bc)
        _saturation, res_dict = self.downsample(df_cell)

        # summary
        self.get_summary(CB_describe, CB_total_Genes,
                         CB_reads_count, reads_mapped_to_transcriptome)

        self.report_prepare()

        self.add_content_item('metric', downsample_summary=res_dict)
        self._clean_up()


@utils.add_log
def count_capture_rna(args):
    # TODO!
    # need barcode_capture_rna
    with Count_capture_rna(args) as runner:
        runner.run()


def get_opts_count_capture_rna(parser, sub_program):
    get_opts_count(parser, sub_program)
