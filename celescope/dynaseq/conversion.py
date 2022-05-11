#!/usr/bin/env python
# v1.0

import pysam
import os
import subprocess
import numpy as np
import pandas as pd
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT

class Conversion(Step):
    """
    ## Features
    - Get conversion pos in each read.
        - Get snp info. 

    ## Output
    - `{sample}.PosTag.bam` Bam file with conversion info.
    - `{sample}.PosTag.csv` SNP info in csv format.
    """

    def __init__(self, args,display_title=None):
        
        super().__init__(args,display_title)

        self.display_title = display_title
        # input files
        self.ifile = os.path.join(args.outdir, args.sample+'.bam')
        self.strandednessfile = args.strand
        self.inbam = args.bam
        self.bcfile = args.cell 
        self.thread = args.thread

        # output files
        self.outfile_bam = os.path.join(args.outdir, args.sample+'.PosTag.bam')
        self.outfile_csv = os.path.join(args.outdir, args.sample+'.PosTag.csv')

    @utils.add_log
    def run(self):
        ##Filter and sort
        self.fltSort(self.inbam, self.ifile, self.bcfile, self.thread)
        cmd = ['samtools index', self.ifile]
        self.run_cmd(cmd)

        # Adding tags
        self.addTags(self.ifile, self.outfile_bam, self.strandednessfile)
        cmd = ['samtools index', self.outfile_bam]
        self.run_cmd(cmd)

        # Obtaining conversion positions
        bam = pysam.AlignmentFile(self.outfile_bam, 'rb')
        ContigLocs, AnnoteLocs = self.CountConvperPos(bam)

        # Obtaining coverage over conversion position
        ConvsPerPos, CoverofPosWithConvs = self.CountReadConverPerConvPos(bam, ContigLocs)
        A = self.ExportasVcf(ConvsPerPos, CoverofPosWithConvs, AnnoteLocs)
        A['sample'] = self.sample
        # Saving result
        A.to_csv(self.outfile_csv)
        bam.close()

        cmd = ['rm', self.ifile]
        self.run_cmd(cmd)
        cmd = ['rm', self.ifile+'.bai']
        self.run_cmd(cmd)

    def run_cmd(self, cmd):
        subprocess.call(' '.join(cmd), shell=True)

    @utils.add_log
    def CountConvperPos(self, bamfile):
        ContigLocs = {}
        AnnoteLocs = {}
        for read in bamfile.fetch():
            try:
                if read.get_tag('ST') == '+':
                    locs = read.get_tag('TL')
                else:
                    locs = read.get_tag('AL')
                if locs[0] != 0:
                    if read.reference_name in ContigLocs:
                        ContigLocs[read.reference_name].extend(locs)
                    else:
                        ContigLocs[read.reference_name] = list(locs)
                    if read.reference_name not in AnnoteLocs:
                        for i, each in enumerate(locs):
                            if i == 0:
                                AnnoteLocs[read.reference_name] = {each: read.get_tag('XT')}
                            else:
                                AnnoteLocs[read.reference_name][each] = read.get_tag('XT')
                    else:
                        for i, each in enumerate(locs):
                            if each not in AnnoteLocs[read.reference_name]:
                                AnnoteLocs[read.reference_name][each] = read.get_tag('XT')
            except (ValueError, KeyError):
                continue
        return ContigLocs, AnnoteLocs

    @utils.add_log
    def CountReadConverPerConvPos(self, bam, ContigLocs):
        ConvsPerPos = {}
        CoverofPosWithConvs = {}
        for key in ContigLocs.keys():
            ContigLocs[key] = sorted(ContigLocs[key])
            ConvsPerPos[key] = {}
            k = 0
            current = ContigLocs[key][k]
            k += 1
            nextone = ContigLocs[key][k]
            while k < len(ContigLocs[key])-1:
                ConvsPerPos[key][current] = 1
                while current == nextone and k < len(ContigLocs[key])-1:
                    k += 1
                    nextone = ContigLocs[key][k]
                    ConvsPerPos[key][current] += 1
                current = nextone
                if k < len(ContigLocs[key])-1:
                    k += 1
                    nextone = ContigLocs[key][k]

            CoverofPosWithConvs[key] = {}
            for key2 in ConvsPerPos[key].keys():
                try:
                    CoverofPosWithConvs[key][key2] = bam.count(key, key2, key2+1)
                except ValueError:
                    continue
        return ConvsPerPos, CoverofPosWithConvs

    @utils.add_log
    def ExportasVcf(self, ConvsPerPos, CoverofPosWithConvs, AnnoteLocs):
        #Chrom, Pos , ConvsPerPs, CoverofPosWithConvs
        Outputdf = pd.DataFrame(columns=['pos2', 'convs', 'covers', 'chrom', 'posratio'])
        for key in ConvsPerPos.keys():
            df = pd.DataFrame.from_dict(ConvsPerPos[key], orient='index')
            df1 = pd.DataFrame.from_dict(CoverofPosWithConvs[key], orient='index')
            df.index.name = 'pos'
            df1.index.name = 'pos'
            df.columns = ['convs']
            df1.columns = ['covers']
            df2 = df.join(df1)
            df2['pos2'] = df2.index
            df2.index = np.arange(df2.shape[0])
            df2['chrom'] = np.repeat(key, df2.shape[0])
            df2['posratio'] = df2['convs']/df2['covers']
            df3 = pd.DataFrame.from_dict(AnnoteLocs[key], orient='index')
            df3.columns = ['gene_id']
            df2 = df2.join(df3, on='pos2')
            Outputdf = Outputdf.append(df2)
        return Outputdf.reset_index(drop=True)

    def createTag(self, d):
        return ''.join([''.join(key) + str(d[key]) + ';' for key in d.keys()])[:-1]

    def convInRead(self, read, qual=20):
        specific_conversions = {}
        total_content = {'a': 0, 'c': 0, 'g': 0, 't': 0}
        specific_conversions[('c', 'A')] = 0
        specific_conversions[('g', 'A')] = 0
        specific_conversions[('t', 'A')] = 0
        specific_conversions[('a', 'C')] = 0
        specific_conversions[('g', 'C')] = 0
        specific_conversions[('t', 'C')] = 0
        specific_conversions[('a', 'G')] = 0
        specific_conversions[('c', 'G')] = 0
        specific_conversions[('t', 'G')] = 0
        specific_conversions[('a', 'T')] = 0
        specific_conversions[('c', 'T')] = 0
        specific_conversions[('g', 'T')] = 0
        specific_conversions[('a', 'N')] = 0
        specific_conversions[('c', 'N')] = 0
        specific_conversions[('g', 'N')] = 0
        specific_conversions[('t', 'N')] = 0

        tC_loc = []
        aG_loc = []

        try:
            refseq = read.get_reference_sequence().lower()
        except (UnicodeDecodeError):
            refseq = ''

        for base in total_content.keys():
            total_content[base] += refseq.count(base)
        for pair in read.get_aligned_pairs(with_seq=True):
            try:
                if pair[0] is not None and pair[1] is not None and pair[2] is not None:
                    if str(pair[2]).islower() and not read.query_qualities[pair[0]] < qual:
                        specific_conversions[(pair[2], read.seq[pair[0]])] += 1
                        if (pair[2], read.seq[pair[0]]) == ('t', 'C'):
                            tC_loc.append(pair[1])
                        if (pair[2], read.seq[pair[0]]) == ('a', 'G'):
                            aG_loc.append(pair[1])
            except (UnicodeDecodeError, KeyError):
                continue
        SC_tag = self.createTag(specific_conversions)
        TC_tag = self.createTag(total_content)

        if len(tC_loc) == 0:
            tC_loc.append(0)
        if len(aG_loc) == 0:
            aG_loc.append(0)
        return SC_tag, TC_tag, tC_loc, aG_loc

    @utils.add_log
    def addTags(self, bamfilename, outputname, strandednessfile):
        bamfile = pysam.AlignmentFile(bamfilename, 'rb')
        mod_bamfile = pysam.AlignmentFile(outputname, mode='wb', template=bamfile)
        strandedness = pd.read_csv(strandednessfile, header=None, index_col=0)
        for read in bamfile.fetch():
            try:
                tags = self.convInRead(read)
                read.set_tag('SC', tags[0], 'Z')
                read.set_tag('TC', tags[1], 'Z')
                read.set_tag('TL', tags[2])
                read.set_tag('AL', tags[3])
                read.set_tag('ST', strandedness.loc[read.get_tag('XT')][1])
                mod_bamfile.write(read)
            except (ValueError, KeyError):
                continue

        bamfile.close()
        mod_bamfile.close()

    @utils.add_log
    def fltSort(self, bamfilename, outfile_bam, cellfile, thread=8):
        bamfile = pysam.AlignmentFile(bamfilename, 'rb')
        mod_bamfile = pysam.AlignmentFile(outfile_bam, mode='wb', template=bamfile)
        cells = {}
        with open(cellfile) as f:
            for i in f:
                cells[i.strip()] = 1
        for read in bamfile.fetch(until_eof=True):
            try:
                if not read.has_tag('GX'):
                    continue
                if read.get_tag("CB") not in cells:
                    continue
                mod_bamfile.write(read)
            except (ValueError, KeyError):
                continue
        bamfile.close()
        mod_bamfile.close()

        cmd = ['samtools sort -@', str(thread), '-o', outfile_bam+'.bam', outfile_bam]
        self.run_cmd(cmd)
        cmd = ['mv', outfile_bam+'.bam', outfile_bam]
        self.run_cmd(cmd)


@utils.add_log
def conversion(args):

    with Conversion(args,display_title='Conversion') as runner:
        runner.run()


def get_opts_conversion(parser, sub_program):
    parser.add_argument('--strand', help=HELP_DICT['strand'], required=True)
    if sub_program:
        parser.add_argument(
            "--bam", help=HELP_DICT['bam_for_conversion'], required=True)
        parser.add_argument("--cell", help=HELP_DICT['cell'], required=True)
        parser = s_common(parser)
    return parser
