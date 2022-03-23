import glob
import gzip
import importlib
import logging
import os
import re
import resource
import subprocess
import time
from collections import Counter, defaultdict
from datetime import timedelta
from functools import wraps

import pandas as pd
import pysam

import celescope.tools
from celescope.tools.__init__ import FILTERED_MATRIX_DIR_SUFFIX, BARCODE_FILE_NAME 
from celescope.__init__ import ROOT_PATH


def add_log(func):
    '''
    logging start and done.
    '''
    logFormatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper


def using(point=""):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
        ''' % (point, usage[0], usage[1],
               usage[2]/1024.0)


def add_mem(func):
    '''
    logging mem.
    '''
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info(using("before"))
        result = func(*args, **kwargs)
        logger.info(using("after"))
        return result

    wrapper.logger = logger
    return wrapper



def generic_open(file_name, *args, **kwargs):
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, *args, **kwargs)
    else:
        file_obj = open(file_name, *args, **kwargs)
    return file_obj


class Gtf_dict(dict):
    '''

    key: gene_id
    value: gene_name
    '''

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

        Returns:
            {gene_id: gene_name} dict
        """

        gene_id_pattern = re.compile(r'gene_id "(\S+)";')
        gene_name_pattern = re.compile(r'gene_name "(\S+)"')
        id_name = {}
        c = Counter()
        with generic_open(self.gtf_file) as f:
            for line in f:
                if not line.strip():
                    continue
                if line.startswith('#'):
                    continue
                tabs = line.split('\t')
                gtf_type, attributes = tabs[2], tabs[-1]
                if gtf_type == 'gene':
                    gene_id = gene_id_pattern.findall(attributes)[-1]
                    gene_names = gene_name_pattern.findall(attributes)
                    if not gene_names:
                        gene_name = gene_id
                    else:
                        gene_name = gene_names[-1]
                    c[gene_name] += 1
                    if c[gene_name] > 1:
                        if gene_id in id_name:
                            assert id_name[gene_id] == gene_name, (
                                'one gene_id with multiple gene_name '
                                f'gene_id: {gene_id}, '
                                f'gene_name this line: {gene_name}'
                                f'gene_name previous line: {id_name[gene_id]}'
                            )
                            self.load_gtf.logger.warning(
                                'duplicated (gene_id, gene_name)'
                                f'gene_id: {gene_id}, '
                                f'gene_name {gene_name}'
                            )
                            c[gene_name] -= 1
                        else:
                            gene_name = f'{gene_name}_{c[gene_name]}'
                    id_name[gene_id] = gene_name
        self.update(id_name)

    def __getitem__(self, key):
        '''if key not exist, return key'''
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


def get_bed_file_path(panel):
    bed_file_path = f'{ROOT_PATH}/data/snp/panel/{panel}.bed'
    if not os.path.exists(bed_file_path):
        return None
    else:
        return bed_file_path


def get_gene_region_from_bed(panel):
    """
    Returns 
    - genes
    - position_df with 'Chromosome', 'Start', 'End'
    """
    file_path = get_bed_file_path(panel)
    bed_file_df = pd.read_table(file_path,
                                usecols=[0, 1, 2, 3],
                                names=['Chromosome', 'Start', 'End', 'Gene'],
                                sep="\t")
    position_df = bed_file_df.loc[:, ['Chromosome', 'Start', 'End']]
    genes = set(bed_file_df.loc[:, 'Gene'].to_list())
    return genes, position_df


def read_fasta(fasta_file, equal=False):
    # seq must have equal length
    fa_dict = {}
    length = None
    with pysam.FastxFile(fasta_file) as infile:
        for index, record in enumerate(infile):
            seq = record.sequence
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
    if (length != length2):
        raise Exception(f"string1({length}) and string2({length2}) do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance


def format_number(number: int) -> str:
    return format(number, ",")


def genDict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: genDict(dim - 1, valType=valType))



def parse_annovar(annovar_file):
    df = pd.DataFrame(columns=['Gene', 'mRNA', 'Protein', 'COSMIC'])
    with open(annovar_file, 'rt') as f:
        index = 0
        for line in f:
            index += 1
            if index == 1:
                continue
            attrs = line.split('\t')
            gene = attrs[6]
            func = attrs[5]
            if func == 'exonic':
                changes = attrs[9]
                cosmic = attrs[10]
            else:
                changes = attrs[7]
                cosmic = attrs[8]
            change_list = list()
            for change in changes.split(','):
                change_attrs = change.split(':')
                mRNA = ''
                protein = ''
                for change_index in range(len(change_attrs)):
                    change_attr = change_attrs[change_index]
                    if change_attr.startswith('c.'):
                        base = change_attr.strip('c.')
                        exon = change_attrs[change_index - 1]
                        mRNA = f'{exon}:{base}'
                    if change_attr.startswith('p.'):
                        protein = change_attr.strip('p.')
                if not (mRNA, protein) in change_list:
                    change_list.append((mRNA, protein))
            combine = [','.join(item) for item in list(zip(*change_list))]
            mRNA = combine[0]
            protein = combine[1]
            df = df.append({
                'Gene': gene,
                'mRNA': mRNA,
                'Protein': protein,
                'COSMIC': cosmic,
            }, ignore_index=True)
    return df

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
        raise TypeError('pattern_list must be a list')

    match_list = []
    for pattern in pattern_list:
        files = glob.glob(pattern)
        if files:
            for f in files:
                match_list.append(f)
    
    if len(match_list) == 0:
        raise FileNotFoundError(f'No file found for {pattern_list}')
    
    if len(match_list) > 1:
        raise MultipleFileFoundError(
            f'More than one file found for pattern: {pattern_list}\n'
            f'File found: {match_list}'
        )
    
    return match_list[0]


@add_log
def get_barcode_from_matrix_dir(matrix_dir):
    """
    Returns:
        match_barcode: list
        no_match_barcode: int
    """
    barcode_file_pattern_list = []
    for barcode_file_name in BARCODE_FILE_NAME:
        barcode_file_pattern_list.append(f"{matrix_dir}/{barcode_file_name}")
  
    match_barcode_file = glob_file(barcode_file_pattern_list)
    get_barcode_from_matrix_dir.logger.info(f"Barcode file:{match_barcode_file}")
    match_barcode, n_match_barcode = read_one_col(match_barcode_file)

    return match_barcode, n_match_barcode


@add_log
def get_matrix_dir_from_match_dir(match_dir):
    """
    Returns:
        matrix_dir: PosixPath object
    """
    matrix_dir_pattern_list = []
    for matrix_dir_suffix in FILTERED_MATRIX_DIR_SUFFIX:
        matrix_dir_pattern_list.append(f"{match_dir}/*count/*{matrix_dir_suffix}")
  
    matrix_dir = glob_file(matrix_dir_pattern_list)
    get_matrix_dir_from_match_dir.logger.info(f"Matrix_dir :{matrix_dir}")

    return matrix_dir


@add_log
def get_barcode_from_match_dir(match_dir):
    '''
    multi version compatible
    Returns:
        match_barcode: list
        no_match_barcode: int
    '''
    matrix_dir = get_matrix_dir_from_match_dir(match_dir)
    return get_barcode_from_matrix_dir(matrix_dir)


@add_log
def parse_match_dir(match_dir):
    '''
    return dict
    keys: 'match_barcode', 'n_match_barcode', 'matrix_dir', 'tsne_coord'
    '''
    match_dict = {}

    pattern_dict = {
        'tsne_coord': [f'{match_dir}/*analysis*/*tsne_coord.tsv'],
        'markers': [f'{match_dir}/*analysis*/*markers.tsv'],
    }

    for file_key in pattern_dict:
        file_pattern= pattern_dict[file_key]
        try:
            match_file = glob_file(file_pattern)
        except FileNotFoundError:
            parse_match_dir.logger.warning(f"No {file_key} found in {match_dir}")
        else:
            match_dict[file_key] = match_file

    match_dict['matrix_dir'] = get_matrix_dir_from_match_dir(match_dir)
    match_barcode, n_match_barcode = get_barcode_from_match_dir(match_dir)
    match_dict['match_barcode'] = match_barcode
    match_dict['n_match_barcode'] = n_match_barcode

    return match_dict




def get_scope_bc(bctype):
    root_path = os.path.dirname(celescope.__file__)
    linker_f = glob.glob(f'{root_path}/data/chemistry/{bctype}/linker*')[0]
    whitelist_f = f'{root_path}/data/chemistry/{bctype}/bclist'
    return linker_f, whitelist_f


def fastq_line(name, seq, qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'


def find_assay_init(assay):
    init_module = importlib.import_module(f"celescope.{assay}.__init__")
    return init_module


def find_step_module(assay, step):
    init_module = find_assay_init(assay)
    try:
        step_module = importlib.import_module(f"celescope.{assay}.{step}")
    except ModuleNotFoundError:
        try:
            step_module = importlib.import_module(f"celescope.tools.{step}")
        except ModuleNotFoundError:
            module_path = init_module.IMPORT_DICT[step]
            step_module = importlib.import_module(f"{module_path}.{step}")

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
            folder = 'tools'
        except ModuleNotFoundError:
            module_path = init_module.IMPORT_DICT[step]
            step_module = importlib.import_module(f"{module_path}.{step}")
            folder = module_path.split('.')[1]

    return step_module, folder


def sort_bam(input_bam, output_bam, threads=1):
    cmd = (
        f'samtools sort {input_bam} '
        f'-o {output_bam} '
        f'--threads {threads} '
    )
    subprocess.check_call(cmd, shell=True)


def index_bam(input_bam):
    cmd = f"samtools index {input_bam}"
    subprocess.check_call(cmd, shell=True)


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")


class Samtools():
    def __init__(self, in_bam, out_bam, threads=1, debug=False):
        self.in_bam = in_bam
        self.out_bam = out_bam
        self.threads = threads
        self.temp_sam_file = f"{self.out_bam}_sam.temp"
        self.debug = debug

    @add_log
    def samtools_sort(self, in_file, out_file, by='coord'):
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

    def sort_bam(self, by='coord'):
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
            with pysam.AlignmentFile(self.temp_sam_file, "w", header=header) as temp_sam:
                for read in original_bam:
                    attr = read.query_name.split('_')
                    barcode = attr[0]
                    umi = attr[1]
                    read.set_tag(tag='CB', value=barcode, value_type='Z')
                    read.set_tag(tag='UB', value=umi, value_type='Z')
                    # assign to some gene
                    if read.has_tag('XT'):
                        gene_id = read.get_tag('XT')
                        # if multi-mapping reads are included in original bam,
                        # there are multiple gene_ids
                        if ',' in gene_id:
                            gene_name = [gtf_dict[i] for i in gene_id.split(',')]
                            gene_name = ','.join(gene_name)
                        else:
                            gene_name = gtf_dict[gene_id]
                        read.set_tag(tag='GN', value=gene_name, value_type='Z')
                        read.set_tag(tag='GX', value=gene_id, value_type='Z')
                    temp_sam.write(read)

    @add_log
    def add_RG(self, barcodes):
        """
        barcodes list
        """

        with pysam.AlignmentFile(self.in_bam, "rb") as original_bam:
            header = original_bam.header.to_dict()
            header['RG'] = []
            for index, barcode in enumerate(barcodes):
                header['RG'].append({
                    'ID': barcode,
                    'SM': index + 1,
                })

            with pysam.AlignmentFile(self.temp_sam_file, "w", header=header) as temp_sam:
                for read in original_bam:
                    read.set_tag(tag='RG', value=read.get_tag('CB'), value_type='Z')
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
    df_index = pd.read_csv(CID_file, sep='\t', index_col=0).reset_index()
    df_valid = df_index[df_index['valid'] == True]
    return df_index, df_valid


def get_assay_text(assay):
    """
    add sinlge cell prefix
    """
    return 'Single-cell ' + assay


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
    if arg_value and arg_value.strip() != 'None':
        return True
    else:
        return False


def parse_vcf_to_df(vcf_file, cols=('chrom', 'pos', 'alleles'), infos=('VID', 'CID')):
    """
    Read cols and infos into pandas df
    """
    vcf = pysam.VariantFile(vcf_file)
    df = pd.DataFrame(columns=[col.capitalize() for col in cols] + infos)
    rec_dict = {}
    for rec in vcf.fetch():

        for col in cols:
            rec_dict[col.capitalize()] = getattr(rec, col)
            if col == 'alleles':
                rec_dict['Alleles'] = '-'.join(rec_dict['Alleles'])

        for info in infos:
            rec_dict[info] = rec.info[info]

        '''
        rec_dict['GT'] = [s['GT'] for s in rec.samples.values()][0]
        rec_dict['GT'] = [str(item) for item in rec_dict['GT']]
        rec_dict['GT'] = '/'.join(rec_dict['GT'])
        '''

        df = df.append(pd.Series(rec_dict), ignore_index=True)
    vcf.close()
    return df