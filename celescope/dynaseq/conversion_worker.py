import subprocess
import pysam
import pandas as pd
from collections import defaultdict


def createTag(d):
    return "".join(["".join(key) + str(d[key]) + ";" for key in d.keys()])[:-1]


def convInRead(read, qual=20):
    specific_conversions = {}
    total_content = {"a": 0, "c": 0, "g": 0, "t": 0}
    specific_conversions[("c", "A")] = 0
    specific_conversions[("g", "A")] = 0
    specific_conversions[("t", "A")] = 0
    specific_conversions[("a", "C")] = 0
    specific_conversions[("g", "C")] = 0
    specific_conversions[("t", "C")] = 0
    specific_conversions[("a", "G")] = 0
    specific_conversions[("c", "G")] = 0
    specific_conversions[("t", "G")] = 0
    specific_conversions[("a", "T")] = 0
    specific_conversions[("c", "T")] = 0
    specific_conversions[("g", "T")] = 0
    specific_conversions[("a", "N")] = 0
    specific_conversions[("c", "N")] = 0
    specific_conversions[("g", "N")] = 0
    specific_conversions[("t", "N")] = 0

    tC_loc = []
    aG_loc = []

    try:
        refseq = read.get_reference_sequence().lower()
    except (UnicodeDecodeError, AssertionError):
        return 0

    for base in total_content.keys():
        total_content[base] += refseq.count(base)
    if read.get_tag("nM") > 0:
        for pair in read.get_aligned_pairs(with_seq=True):
            try:
                if pair[0] is not None and pair[1] is not None and pair[2] is not None:
                    if (
                        str(pair[2]).islower()
                        and not read.query_qualities[pair[0]] < qual
                    ):
                        specific_conversions[(pair[2], read.seq[pair[0]])] += 1
                        if (pair[2], read.seq[pair[0]]) == ("t", "C"):
                            tC_loc.append(pair[1])
                        if (pair[2], read.seq[pair[0]]) == ("a", "G"):
                            aG_loc.append(pair[1])
            except (UnicodeDecodeError, KeyError):
                continue

    SC_tag = createTag(specific_conversions)
    TC_tag = createTag(total_content)

    if len(tC_loc) == 0:
        tC_loc.append(0)
    if len(aG_loc) == 0:
        aG_loc.append(0)
    return SC_tag, TC_tag, tC_loc, aG_loc


def process_bam(bamfilename, tmpoutbam, cells, strandedness, qual):
    site_depth = defaultdict(int)

    save = pysam.set_verbosity(0)
    bamfile = pysam.AlignmentFile(bamfilename, "rb")
    header = bamfile.header
    mod_bamfile = pysam.AlignmentFile(
        tmpoutbam, mode="wb", header=header, check_sq=False
    )
    pysam.set_verbosity(save)

    class GeneError(Exception):
        pass

    for read in bamfile.fetch(until_eof=True):
        try:
            if (not read.has_tag("GX")) or read.get_tag("GX") == "-":
                continue
            if read.get_tag("GX") not in strandedness:
                raise GeneError
            if read.get_tag("CB") not in cells:
                continue

            tags = convInRead(read, qual)
            if tags == 0:
                continue
            read.set_tag("SC", tags[0], "Z")
            read.set_tag("TC", tags[1], "Z")
            read.set_tag("TL", tags[2])
            read.set_tag("AL", tags[3])
            read.set_tag("ST", strandedness[read.get_tag("GX")])
            mod_bamfile.write(read)

            if strandedness[read.get_tag("GX")] == "+":
                locs = tags[2]
            else:
                locs = tags[3]
            if locs[0] != 0:
                for loc in locs:
                    site = f"{read.reference_name}+{loc}"
                    site_depth[site] += 1

        except (ValueError, KeyError):
            continue
        except GeneError:
            print(f"{read.get_tag('GX')} is not in gtf file, please check your files.")
            continue

    bamfile.close()
    mod_bamfile.close()

    if len(site_depth) == 0:  ## if no conversion site detected
        return pd.DataFrame()

    df = pd.DataFrame.from_dict(site_depth, orient="index")
    df.columns = ["convs"]
    return df


def count_read_cover_per_conv_pos(outbam, df_conv):
    bamfile = outbam
    df = df_conv
    cover_of_pos_with_convs = {}
    if df.shape[0] == 0:
        return cover_of_pos_with_convs
    df = df.reset_index()
    df[["chrom", "pos"]] = df["index"].str.split("+", expand=True)
    df["pos"] = df["pos"].astype(int)
    try:
        cmd = f"samtools index {bamfile} 2>&1 "
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        cmd = f"samtools index -c {bamfile} 2>&1 "
        subprocess.check_call(cmd, shell=True)
    save = pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(bamfile, "rb")
    pysam.set_verbosity(save)
    contig_locs = df.groupby("chrom")["pos"].apply(list).to_dict()
    for key in contig_locs.keys():
        contig_locs[key] = sorted(contig_locs[key])
        cover_of_pos_with_convs[key] = {}
        for key2 in contig_locs[key]:
            try:
                cover_of_pos_with_convs[key][key2] = bam.count(key, key2, key2 + 1)
            except ValueError:
                continue
    bam.close()
    return cover_of_pos_with_convs


def conv_candidate(df_conv, df_cover):
    df = df_conv
    if df.shape[0] == 0:
        df = pd.DataFrame(columns=["chrom", "pos", "convs", "covers"])
        df = df.set_index(["chrom", "pos"])
    else:
        df = df.reset_index()
        df[["chrom", "pos"]] = df["index"].str.split("+", expand=True)
        df["pos"] = df["pos"].astype(int)
        df = df.set_index(["chrom", "pos"])
        df.drop(["index"], axis=1, inplace=True)
        dep = pd.DataFrame.from_dict(df_cover, orient="index")
        dep = dep.reset_index()
        dep = dep.melt(id_vars="index").dropna()
        dep.columns = ["chrom", "pos", "covers"]
        dep = dep.set_index(["chrom", "pos"])
        df = pd.concat([df, dep], axis=1)

    return df


def conversion_process(args):
    bamfilename, tmpoutbam, cells, strandedness, qual = args
    df_conv = process_bam(bamfilename, tmpoutbam, cells, strandedness, qual)
    df_cover = count_read_cover_per_conv_pos(tmpoutbam, df_conv)
    df = conv_candidate(df_conv, df_cover)
    return df
