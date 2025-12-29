import argparse
import glob
import gzip
import importlib
import json
import logging
import os
import re
import subprocess
import sys
import time
import unittest
from collections import Counter, OrderedDict, defaultdict
from datetime import timedelta
from functools import wraps
from typing import List

import pandas as pd
import pysam

from celescope.__init__ import ROOT_PATH
from celescope.tools.__init__ import (
    BARCODE_FILE_NAME,
    FILTERED_MATRIX_DIR_SUFFIX,
    OUTS_DIR,
)


def add_log(func):
    """
    logging start and done.
    """
    logFormatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    module = func.__module__
    name = func.__name__
    logger_name = f"{module}.{name}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler(sys.stderr)
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info("start...")
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info("done. time used: %s", used)
        return result

    wrapper.logger = logger
    return wrapper


def generic_open(file_name, *args, **kwargs):
    if file_name.endswith(".gz"):
        file_obj = gzip.open(file_name, *args, **kwargs)
    else:
        file_obj = open(file_name, *args, **kwargs)
    return file_obj


class Gtf_dict(dict):
    """
    key: gene_id
    value: gene_name
    If the key does not exist, return key. This is to avoid the error:
        The gtf file contains one exon lines with a gene_id, but do not contain a gene line with the same gene_id. FeatureCounts
        work correctly under this condition, but the gene_id will not appear in the Gtf_dict.
    """

    def __init__(self, gtf_file):
        super().__init__()
        self.gtf_file = gtf_file
        self.load_gtf()

    @add_log
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
        id_name = {}
        c = Counter()
        with generic_open(self.gtf_file, mode="rt") as f:
            for line in f:
                if not line.strip():
                    continue
                if line.startswith("#"):
                    continue
                tabs = line.split("\t")
                gtf_type, attributes = tabs[2], tabs[-1]
                if gtf_type == "gene":
                    try:
                        gene_id = gene_id_pattern.findall(attributes)[-1]
                    except IndexError:
                        print(line)
                    gene_names = gene_name_pattern.findall(attributes)
                    if not gene_names:
                        gene_name = gene_id
                    else:
                        gene_name = gene_names[-1]
                    c[gene_name] += 1
                    if c[gene_name] > 1:
                        if gene_id in id_name:
                            assert id_name[gene_id] == gene_name, (
                                "one gene_id with multiple gene_name "
                                f"gene_id: {gene_id}, "
                                f"gene_name this line: {gene_name}"
                                f"gene_name previous line: {id_name[gene_id]}"
                            )
                            self.load_gtf.logger.warning(
                                "duplicated (gene_id, gene_name)"
                                f"gene_id: {gene_id}, "
                                f"gene_name {gene_name}"
                            )
                            c[gene_name] -= 1
                        else:
                            gene_name = f"{gene_name}_{c[gene_name]}"
                    id_name[gene_id] = gene_name
        self.update(id_name)

    def __getitem__(self, key):
        """if key not exist, return key"""
        return dict.get(self, key, key)


def read_one_col(file):
    """
    Read file with one column. Strip each line.
    Returns col_list, line number
    """
    df = pd.read_csv(file, header=None)
    col1 = list(df.iloc[:, 0])
    col1 = [item.strip() for item in col1]
    num = len(col1)
    return col1, num


def one_col_to_list(file) -> list:
    """
    Read file with one column. Strip each line.
    Returns col_list
    """
    df = pd.read_csv(file, header=None)
    col1 = list(df.iloc[:, 0])
    return [item.strip() for item in col1]


def two_col_to_dict(file):
    """
    Read file with two columns.
    Returns dict
    """
    df = pd.read_csv(file, header=None, sep="\t")
    df = df.dropna()
    return OrderedDict(zip(df[0], df[1]))


def get_bed_file_path(panel):
    bed_file_path = f"{ROOT_PATH}/data/snp/panel/{panel}.bed"
    if not os.path.exists(bed_file_path):
        return None
    else:
        return bed_file_path


def get_gene_region_from_bed(bed) -> tuple[set, pd.DataFrame]:
    """
    Returns
    - genes: set
    - position_df with 'Chromosome', 'Start', 'End'
    """
    bed_file_df = pd.read_table(
        bed,
        usecols=[0, 1, 2, 3],
        names=["Chromosome", "Start", "End", "Gene"],
        sep="\t",
    )
    position_df = bed_file_df.loc[:, ["Chromosome", "Start", "End"]]
    genes = set(bed_file_df.loc[:, "Gene"].to_list())
    return genes, position_df


def read_fasta(fasta_file, equal=False):
    """
    Args:
        equal: if True, seq in fasta must have equal length
    Returns:
        {seq_id: seq} dict, length
    """
    fa_dict = {}
    length = None
    with pysam.FastxFile(fasta_file) as infile:
        for index, record in enumerate(infile):
            seq = str(record.sequence).upper()
            if index == 0:
                length = len(seq)
            if equal:
                if length != len(seq):
                    raise Exception(f"{fasta_file} have different seq length")
            fa_dict[record.name] = seq
    return fa_dict, length


def hamming_correct(string1, string2):
    threshold = len(string1) / 10 + 1
    if hamming_distance(string1, string2) < threshold:
        return True
    return False


def hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if length != length2:
        raise Exception(
            f"string1({length}) and string2({length2}) do not have same length"
        )
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance


def format_number(number: int) -> str:
    return format(number, ",")


def nested_defaultdict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: nested_defaultdict(dim - 1, valType=valType))


class MultipleFileFoundError(Exception):
    pass


def glob_file(pattern_list: list):
    """
    glob file among pattern list
    Returns:
        PosixPath object
    Raises:
        FileNotFoundError: if no file found
        MultipleFileFound: if more than one file is found
    """
    if not isinstance(pattern_list, list):
        raise TypeError("pattern_list must be a list")

    match_list = []
    for pattern in pattern_list:
        files = glob.glob(pattern)
        if files:
            for f in files:
                match_list.append(f)

    if len(match_list) == 0:
        raise FileNotFoundError(f"No file found for {pattern_list}")

    if len(match_list) > 1:
        raise MultipleFileFoundError(
            f"More than one file found for pattern: {pattern_list}\n"
            f"File found: {match_list}"
        )

    return match_list[0]


def get_matrix_file_path(matrix_dir, file_name):
    """
    compatible with non-gzip file
    """
    non_gzip = file_name.strip(".gz")
    file_path_list = [f"{matrix_dir}/{file_name}", f"{matrix_dir}/{non_gzip}"]
    for file_path in file_path_list:
        if os.path.exists(file_path):
            return file_path


@add_log
def get_barcode_from_matrix_dir(matrix_dir):
    """
    Returns:
        match_barcode: list
        n_match_barcode: int
    """

    match_barcode_file = get_matrix_file_path(matrix_dir, BARCODE_FILE_NAME)
    match_barcode, n_match_barcode = read_one_col(match_barcode_file)

    return match_barcode, n_match_barcode


@add_log
def get_matrix_dir_from_match_dir(match_dir):
    """
    Returns:
        matrix_dir: PosixPath object
    """
    matrix_dir = f"{match_dir}/{OUTS_DIR}/{FILTERED_MATRIX_DIR_SUFFIX}"
    if not os.path.exists(matrix_dir):
        raise FileNotFoundError(f"{matrix_dir} not found")

    return matrix_dir


@add_log
def get_barcode_from_match_dir(match_dir):
    """
    multi version compatible
    Returns:
        match_barcode: list
        no_match_barcode: int
    """
    matrix_dir = get_matrix_dir_from_match_dir(match_dir)
    return get_barcode_from_matrix_dir(matrix_dir)


@add_log
def parse_match_dir(match_dir):
    """
    return dict
    keys: 'match_barcode', 'n_match_barcode', 'matrix_dir', 'tsne_coord'
    """
    match_dict = {}

    pattern_dict = {
        "tsne_coord": [f"{match_dir}/{OUTS_DIR}/tsne_coord.tsv"],
        "markers": [f"{match_dir}/{OUTS_DIR}/markers.tsv"],
        "h5ad": [f"{match_dir}/{OUTS_DIR}/rna.h5ad"],
    }

    for file_key in pattern_dict:
        file_pattern = pattern_dict[file_key]
        try:
            match_file = glob_file(file_pattern)
        except FileNotFoundError:
            parse_match_dir.logger.warning(f"No {file_key} found in {match_dir}")
        else:
            match_dict[file_key] = match_file

    match_dict["matrix_dir"] = get_matrix_dir_from_match_dir(match_dir)
    match_barcode, n_match_barcode = get_barcode_from_match_dir(match_dir)
    match_dict["match_barcode"] = match_barcode
    match_dict["n_match_barcode"] = n_match_barcode
    return match_dict


def read_tsne(tsne_file):
    df = pd.read_csv(tsne_file, sep="\t", index_col=0)
    df.index.names = ["barcode"]
    return df


def fastq_line(name, seq, qual):
    return f"@{name}\n{seq}\n+\n{qual}\n"


def fasta_line(name, seq):
    return f">{name}\n{seq}\n"


def find_assay_init(assay):
    init_module = importlib.import_module(f"celescope.{assay}.__init__")
    return init_module


def find_step_module(assay, step):
    file_path_dict = {
        "assay": f"{ROOT_PATH}/{assay}/{step}.py",
        "tools": f"{ROOT_PATH}/tools/{step}.py",
    }

    init_module = find_assay_init(assay)
    if os.path.exists(file_path_dict["assay"]):
        step_module = importlib.import_module(f"celescope.{assay}.{step}")
    elif hasattr(init_module, "IMPORT_DICT") and step in init_module.IMPORT_DICT:
        module_path = init_module.IMPORT_DICT[step]
        step_module = importlib.import_module(f"{module_path}.{step}")
    elif os.path.exists(file_path_dict["tools"]):
        step_module = importlib.import_module(f"celescope.tools.{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {assay}.{step}")

    return step_module


def find_step_module_with_folder(assay, step):
    init_module = find_assay_init(assay)
    folder = ""
    try:
        step_module = importlib.import_module(f"celescope.{assay}.{step}")
        folder = assay
    except ModuleNotFoundError:
        try:
            step_module = importlib.import_module(f"celescope.tools.{step}")
            folder = "tools"
        except ModuleNotFoundError:
            module_path = init_module.IMPORT_DICT[step]
            step_module = importlib.import_module(f"{module_path}.{step}")
            folder = module_path.split(".")[1]

    return step_module, folder


def sort_bam(input_bam, output_bam, threads=1, by="pos"):
    cmd = (
        f"samtools sort {input_bam} "
        f"-o {output_bam} "
        f"--threads {threads} "
        "2>&1 "
    )
    if by == "name":
        cmd += " -n"
    subprocess.check_call(cmd, shell=True)


def index_bam(input_bam):
    cmd = f"samtools index {input_bam} 2>&1 "
    subprocess.check_call(cmd, shell=True)


def add_tag(seg, id_name, correct_dict):
    """
    Args:
        seg: pysam bam segment
        id_name: {gene_id: gene_name}
        correct_dict: {low_seq: high_seq}

    Returns:
        seg with tag added

    """
    attr = seg.query_name.split(":")
    barcode = attr[0]
    umi = attr[1]
    seg.set_tag(tag="CB", value=barcode, value_type="Z")
    if umi in correct_dict:
        umi = correct_dict[umi]
    seg.set_tag(tag="UB", value=umi, value_type="Z")
    # assign to some gene
    if seg.has_tag("XT"):
        gene_id = seg.get_tag("XT")
        # if multi-mapping reads are included in original bam,
        # there are multiple gene_ids
        if "," in gene_id:
            gene_name = [id_name[i] for i in gene_id.split(",")]
            gene_name = ",".join(gene_name)
        else:
            gene_name = id_name[gene_id]
        seg.set_tag(tag="GN", value=gene_name, value_type="Z")
        seg.set_tag(tag="GX", value=gene_id, value_type="Z")

    return seg


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")


class Samtools:
    def __init__(self, in_bam, out_bam, threads=1, debug=False):
        self.in_bam = in_bam
        self.out_bam = out_bam
        self.threads = threads
        self.temp_sam_file = f"{self.out_bam}_sam.temp"
        self.debug = debug

    @add_log
    def samtools_sort(self, in_file, out_file, by="pos"):
        cmd = f"samtools sort {in_file} -o {out_file} --threads {self.threads}"
        if by == "name":
            cmd += " -n"
        self.samtools_sort.logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)

    @add_log
    def samtools_index(self, in_file):
        cmd = f"samtools index {in_file}"
        self.samtools_index.logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)

    def sort_bam(self, by="coord"):
        """sort in_bam"""
        self.samtools_sort(self.in_bam, self.out_bam, by=by)

    def index_bam(self):
        """index out_bam"""
        self.samtools_index(self.out_bam)

    @add_log
    def add_tag(self, gtf_file):
        """
        - CB cell barcode
        - UB UMI
        - GN gene name
        - GX gene id
        """
        gtf_dict = Gtf_dict(gtf_file)

        with pysam.AlignmentFile(self.in_bam, "rb") as original_bam:
            header = original_bam.header
            with pysam.AlignmentFile(
                self.temp_sam_file, "w", header=header
            ) as temp_sam:
                for read in original_bam:
                    attr = read.query_name.split(":")
                    barcode = attr[0]
                    umi = attr[1]
                    read.set_tag(tag="CB", value=barcode, value_type="Z")
                    read.set_tag(tag="UB", value=umi, value_type="Z")
                    # assign to some gene
                    if read.has_tag("XT"):
                        gene_id = read.get_tag("XT")
                        # if multi-mapping reads are included in original bam,
                        # there are multiple gene_ids
                        if "," in gene_id:
                            gene_name = [gtf_dict[i] for i in gene_id.split(",")]
                            gene_name = ",".join(gene_name)
                        else:
                            gene_name = gtf_dict[gene_id]
                        read.set_tag(tag="GN", value=gene_name, value_type="Z")
                        read.set_tag(tag="GX", value=gene_id, value_type="Z")
                    temp_sam.write(read)

    @add_log
    def add_RG(self, barcodes):
        """
        barcodes list
        """

        with pysam.AlignmentFile(self.in_bam, "rb") as original_bam:
            header = original_bam.header.to_dict()
            header["RG"] = []
            for index, barcode in enumerate(barcodes):
                header["RG"].append(
                    {
                        "ID": barcode,
                        "SM": index + 1,
                    }
                )

            with pysam.AlignmentFile(
                self.temp_sam_file, "w", header=header
            ) as temp_sam:
                for read in original_bam:
                    read.set_tag(tag="RG", value=read.get_tag("CB"), value_type="Z")
                    temp_sam.write(read)

    def temp_sam2bam(self, by=None):
        self.samtools_sort(self.temp_sam_file, self.out_bam, by=by)
        self.rm_temp_sam()

    def rm_temp_sam(self):
        cmd = f"rm {self.temp_sam_file}"
        subprocess.check_call(cmd, shell=True)


def read_CID(CID_file):
    """
    return df_index, df_valid
    """
    df_index = pd.read_csv(CID_file, sep="\t", index_col=0).reset_index()
    df_valid = df_index[df_index["valid"] is True]
    return df_index, df_valid


def get_assay_text(assay):
    """
    Deprecated
    add sinlge cell prefix
    deprecated
    """
    return "Single-cell " + assay


def check_arg_not_none(args, arg_name):
    """
    check if args.arg_name is not None
    Args:
        args: argparser args
        arg_name: argparser arg name
    Return:
        bool
    """
    arg_value = getattr(args, arg_name, None)
    if arg_value and arg_value.strip() != "None":
        return True
    else:
        return False


def reverse_complement(dna: str) -> str:
    """Returns the reverse complement of a DNA sequence, allowing 'N' bases.
    >>> reverse_complement("ATCGNTA")
    'TANCGAT'
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(dna))


def get_fastx_read_number(fastx_file):
    """
    get read number using pysam
    """
    n = 0
    with pysam.FastxFile(fastx_file) as f:
        for _ in f:
            n += 1
    return n


@add_log
def dump_dict_to_json(d, json_file):
    with open(json_file, "w") as f:
        json.dump(d, f, indent=4)


@add_log
def barcode_list_stamp(barcode_list, cut=500):
    bc_list, bc_num = read_one_col(barcode_list)
    n, m = 0, 0
    stamp = defaultdict(list)
    for i in bc_list:
        m += 1
        if m <= cut:
            stamp[n].append(i)
        else:
            n += 1
            m = 1
            stamp[n].append(i)
    return stamp, bc_num


def merge_table_files(
    table_files: List[str],
    output_file: str,
    sep: str = "\t",
) -> None:
    if not table_files:
        raise ValueError("Error: empty table_files")

    df = pd.read_csv(table_files[0], sep=sep)

    for file in table_files[1:]:
        temp_df = pd.read_csv(file, sep=sep)
        df = pd.concat([df, temp_df], axis=0)

    df.to_csv(output_file, sep=sep, index=False)


def get_action(parser: argparse.ArgumentParser, option_string: str) -> argparse.Action:
    """
    根据 option_string 找到对应的 argparse.Action 对象。

    参数:
        parser: argparse.ArgumentParser 实例
        option_string: 参数名，例如 '--foo' 或 '-f'

    返回:
        对应的 argparse.Action 对象，如果没找到返回 None
    """
    for action in parser._actions:
        if option_string in action.option_strings:
            return action
    raise ValueError(f"{option_string} not found in parser")


def add_quotes_if_needed(s: str) -> str:
    s = s.strip()
    matches = [" ", "-", ">", "<"]
    if any(char in s for char in matches):
        return f'"{s}"'
    return s


class Test_utils(unittest.TestCase):
    def test_gtf_dict(self):
        import tempfile

        fp = tempfile.NamedTemporaryFile(suffix=".gtf")
        fp.write(
            b'1\tprocessed_transcript\tgene\t11869\t14409\t.\t+\t.\tgene_id "gene_id_with_gene_name"; gene_name "gene_name";\n'
        )
        fp.write(
            b'1\tprocessed_transcript\tgene\t11869\t14409\t.\t+\t.\tgene_id "gene_id_without_gene_name"; \n'
        )
        fp.seek(0)
        gtf_dict = Gtf_dict(fp.name)
        self.assertEqual(gtf_dict["gene_id_with_gene_name"], "gene_name")
        self.assertEqual(
            gtf_dict["gene_id_without_gene_name"], "gene_id_without_gene_name"
        )
        self.assertEqual(gtf_dict["gene_id_not_exist"], "gene_id_not_exist")
        fp.close()


if __name__ == "__main__":
    unittest.main()
