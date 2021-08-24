import argparse
import glob
import gzip
import importlib
import itertools
import json
import logging
import os
import re
import resource
import subprocess
import time
from collections import Counter, defaultdict
from datetime import timedelta
from functools import wraps

import numpy as np
import pandas as pd
import pysam
import xopen

import celescope.tools
from celescope.tools.__init__ import __PATTERN_DICT__
from celescope.capture_virus.otsu import array2hist, makePlot, threshold_otsu

tools_dir = os.path.dirname(celescope.tools.__file__)


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
        if args and hasattr(args[0], 'debug') and args[0].debug:
            logger.setLevel(10)  # debug

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


def arg_str(arg, arg_name):
    '''
    return action store_true arguments as string
    '''
    if arg:
        return '--' + arg_name
    return ''


def format_stat(count, total_count):
    percent = round(count / total_count * 100, 2)
    string = f'{format_number(count)}({percent}%)'
    return string


def multi_opts(assay):
    readme = f'{assay} multi-samples'
    parser = argparse.ArgumentParser(readme)
    parser.add_argument('--mod', help='mod, sjm or shell', choices=['sjm', 'shell'], default='sjm')
    parser.add_argument(
        '--mapfile',
        help='''
            tsv file, 4 columns:
            1st col: LibName;
            2nd col: DataDir;
            3rd col: SampleName;
            4th col: Cell number or match_dir, optional;
        ''',
        required=True)
    parser.add_argument('--chemistry', choices=list(__PATTERN_DICT__.keys()), help='chemistry version', default='auto')
    parser.add_argument('--whitelist', help='cellbarcode list')
    parser.add_argument('--linker', help='linker')
    parser.add_argument('--pattern', help='read1 pattern')
    parser.add_argument('--outdir', help='output dir', default="./")
    parser.add_argument(
        '--adapt',
        action='append',
        help='adapter sequence',
        default=[
            'polyT=A{15}',
            'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'])
    parser.add_argument(
        '--minimum_length',
        dest='minimum_length',
        help='minimum_length',
        default=20)
    parser.add_argument(
        '--nextseq-trim',
        dest='nextseq_trim',
        help='nextseq_trim',
        default=20)
    parser.add_argument(
        '--overlap',
        help='minimum overlap length, default=5',
        default=5)
    parser.add_argument(
        '--lowQual',
        type=int,
        help='max phred of base as lowQual',
        default=0)
    parser.add_argument(
        '--lowNum',
        type=int,
        help='max number with lowQual allowed',
        default=2)
    parser.add_argument(
        '--rm_files',
        action='store_true',
        help='remove redundant fq.gz and bam after running')
    parser.add_argument('--steps_run', help='steps to run', default='all')
    parser.add_argument('--debug', help='debug or not', action='store_true')
    parser.add_argument('--outFilterMatchNmin', help='STAR outFilterMatchNmin', default=0)
    return parser


def link_data(outdir, fq_dict):
    raw_dir = f'{outdir}/data_give/rawdata'
    os.system('mkdir -p %s' % (raw_dir))
    with open(raw_dir + '/ln.sh', 'w') as fh:
        fh.write('cd %s\n' % (raw_dir))
        for s, arr in fq_dict.items():
            fh.write('ln -sf %s %s\n' % (arr[0], s + '_1.fq.gz'))
            fh.write('ln -sf %s %s\n' % (arr[1], s + '_2.fq.gz'))


def generic_open(file_name, *args, **kwargs):
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, *args, **kwargs)
    else:
        file_obj = open(file_name, *args, **kwargs)
    return file_obj


@add_log
def get_id_name_dict(gtf_file):
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
    with generic_open(gtf_file) as f:
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
                        get_id_name_dict.logger.warning(
                            'duplicated (gene_id, gene_name)'
                            f'gene_id: {gene_id}, '
                            f'gene_name {gene_name}'
                        )
                        c[gene_name] -= 1
                    else:
                        gene_name = f'{gene_name}_{c[gene_name]}'
                id_name[gene_id] = gene_name
    return id_name


@add_log
def process_read(
        read2_file, pattern_dict, barcode_dict, linker_dict,
        barcode_length, linker_length):

    # if valid, return (True)
    metrics = defaultdict(int)
    res_dict = genDict(dim=3)
    read2 = xopen.xopen(read2_file, "rt")
    while True:
        line1 = read2.readline()
        line2 = read2.readline()
        read2.readline()
        line4 = read2.readline()
        if not line4:
            break
        metrics['Total Reads'] += 1
        attr = str(line1).strip("@").split("_")
        barcode = str(attr[0])
        umi = str(attr[1])
        seq = line2.strip()
        if linker_length != 0:
            seq_linker = seq_ranges(seq, pattern_dict['L'])
            if len(seq_linker) < linker_length:
                metrics['Reads Unmapped too Short'] += 1
                continue
        if barcode_dict:
            seq_barcode = seq_ranges(seq, pattern_dict['C'])
            if barcode_length != len(seq_barcode):
                miss_length = barcode_length - len(seq_barcode)
                if miss_length > 2:
                    metrics['Reads Unmapped too Short'] += 1
                    continue
                seq_barcode = seq_barcode + "A" * miss_length

        # check linker
        if linker_length != 0:
            valid_linker = False
            for linker_name in linker_dict:
                if hamming_correct(linker_dict[linker_name], seq_linker):
                    valid_linker = True
                    break
        else:
            valid_linker = True

        if not valid_linker:
            metrics['Reads Unmapped Invalid Linker'] += 1
            continue

        # check barcode
        valid_barcode = False
        for barcode_name in barcode_dict:
            if hamming_correct(barcode_dict[barcode_name], seq_barcode):
                res_dict[barcode][barcode_name][umi] += 1
                valid_barcode = True
                break

        if not valid_barcode:
            metrics['Reads Unmapped Invalid Barcode'] += 1
            continue

        # mapped
        metrics['Reads Mapped'] += 1
        if metrics['Reads Mapped'] % 1000000 == 0:
            process_read.logger.info(str(metrics['Reads Mapped']) + " reads done.")

    return res_dict, metrics


def seq_ranges_exception(seq, pattern_dict):
    # get subseq with intervals in arr and concatenate
    length = len(seq)
    for x in pattern_dict:
        if length < x[1]:
            raise Exception(f'invalid seq range {x[0]}:{x[1]} in read')
    return ''.join([seq[x[0]:x[1]]for x in pattern_dict])


def seq_ranges(seq, pattern_dict):
    # get subseq with intervals in arr and concatenate
    return ''.join([seq[x[0]:x[1]]for x in pattern_dict])


def read_one_col(file):
    df = pd.read_csv(file, header=None)
    col1 = list(df.iloc[:, 0])
    num = len(col1)
    return col1, num


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


def gen_stat(df, stat_file):
    # 3cols: item count total_count

    def add_percent(row):
        count = row['count']
        percent = count / row['total_count']
        value = f'{format_number(count)}({round(percent * 100, 2)}%)'
        return value

    df.loc[:, 'value'] = df.loc[:, 'count']
    df.loc[~df['total_count'].isna(), 'value'] = df.loc[~df['total_count'].isna(), :].apply(
        add_percent, axis=1
    )
    df.loc[df['total_count'].isna(), 'value'] = df.loc[df['total_count'].isna(), :].apply(
        lambda row: f'{format_number(row["count"])}', axis=1
    )
    df = df.loc[:, ["item", "value"]]
    df.to_csv(stat_file, sep=":", header=None, index=False)


def get_read(library_id, library_path, read='1'):
    read1_list = [f'_{read}', f'R{read}', f'R{read}_001']
    fq_list = ['fq', 'fastq']
    suffix_list = ["", ".gz"]
    read_pattern_list = [
        f'{library_path}/*{library_id}*{read}.{fq_str}{suffix}'
        for read in read1_list
        for fq_str in fq_list
        for suffix in suffix_list
    ]
    fq_list = [glob.glob(read1_pattern) for read1_pattern in read_pattern_list]
    fq_list = sorted(non_empty for non_empty in fq_list if non_empty)
    fq_list = list(itertools.chain(*fq_list))
    if len(fq_list) == 0:
        print("Allowed R1 patterns:")
        for pattern in read_pattern_list:
            print(pattern)
        raise Exception(
            '\n'
            f'Invalid Read{read} path! \n'
            f'library_id: {library_id}\n'
            f'library_path: {library_path}\n'
        )
    return fq_list


def get_fq(library_id, library_path):
    fq1_list = get_read(library_id, library_path, read='1')
    fq2_list = get_read(library_id, library_path, read='2')
    if len(fq1_list) != len(fq2_list):
        raise Exception("Read1 and Read2 fastq number do not match!")
    fq1 = ",".join(fq1_list)
    fq2 = ",".join(fq2_list)
    return fq1, fq2


def format_number(number: int) -> str:
    return format(number, ",")


def format_metrics(metrics: dict):
    for key in metrics:
        value = metrics[key]
        metrics[key] = format_number(value)
        if isinstance(value, float):
            metrics[key] = round(value, 2)


def format_ratios(ratios: dict):
    for key in ratios:
        ratios[key] = round(ratios[key] * 100, 2)


def get_slope(x, y, window=200, step=10):
    assert len(x) == len(y)
    start = 0
    last = len(x)
    res = [[], []]
    while True:
        end = start + window
        if end > last:
            break
        z = np.polyfit(x[start:end], y[start:end], 1)
        res[0].append(x[start])
        res[1].append(z[0])
        start += step
    return res


def genDict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: genDict(dim - 1, valType=valType))


def report_prepare(outdir, **kwargs):
    json_file = outdir + '/../.data.json'
    if not os.path.exists(json_file):
        data = {}
    else:
        fh = open(json_file)
        data = json.load(fh)
        fh.close()

    for key in kwargs:
        data[key] = kwargs[key]

    with open(json_file, 'w') as fh:
        json.dump(data, fh)


def parse_vcf(vcf_file, cols=('chrom', 'pos', 'alleles'), infos=('VID', 'CID')):
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


def read_barcode_file(match_dir, return_file=False):
    '''
    multi version compatible
    '''
    match_barcode_file1 = glob.glob(
        f"{match_dir}/*count*/*_cellbarcode.tsv")
    match_barcode_file2 = glob.glob(
        f"{match_dir}/*count*/*matrix_10X/*_cellbarcode.tsv")
    match_barcode_file3 = glob.glob(
        f"{match_dir}/*count*/*matrix_10X/*barcodes.tsv")
    match_barcode_file = (
        match_barcode_file1 +
        match_barcode_file2 +
        match_barcode_file3)[0]
    match_barcode, cell_total = read_one_col(match_barcode_file)
    if return_file:
        return match_barcode, (cell_total, match_barcode_file)
    return match_barcode, cell_total


def get_barcodes_from_matrix_dir(matrix_dir):
    barcodes_file = f'{matrix_dir}/barcodes.tsv'
    match_barcode, _cell_total = read_one_col(barcodes_file)
    return match_barcode


def parse_match_dir(match_dir):
    match_dict = {}
    match_barcode, cell_total = read_barcode_file(match_dir)
    match_dict['match_barcode'] = match_barcode
    match_dict['cell_total'] = cell_total
    match_dict['matrix_dir'] = glob.glob(f'{match_dir}/*count*/*matrix_10X')[0]
    match_dict['tsne_coord'] = glob.glob(f'{match_dir}/*analysis*/*tsne_coord.tsv')[0]
    match_dict['markers'] = glob.glob(f'{match_dir}/*analysis*/*markers.tsv')[0]
    try:
        match_dict['rds'] = glob.glob(f'{match_dir}/*analysis/*.rds')[0]
    except IndexError:
        match_dict['rds'] = None
    return match_dict


@add_log
def STAR_util(
    sample,
    outdir,
    input_read,
    genomeDir,
    runThreadN,
    outFilterMatchNmin=35,
    out_unmapped=False,
    outFilterMultimapNmax=1,
    outBAMsortingBinsN=50,
):

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    out_prefix = f'{outdir}/{sample}_'
    out_BAM = out_prefix + "Aligned.sortedByCoord.out.bam"

    cmd = f"STAR \
 --genomeDir {genomeDir} \
 --readFilesIn {input_read}\
 --readFilesCommand zcat\
 --outSAMtype BAM SortedByCoordinate\
 --runThreadN {runThreadN}\
 --outFilterMatchNmin {outFilterMatchNmin}\
 --outFileNamePrefix {out_prefix}\
 --outFilterMultimapNmax {outFilterMultimapNmax}\
 --outBAMsortingBinsN {outBAMsortingBinsN}\
 --limitBAMsortRAM 100000000000 "

    if out_unmapped:
        cmd += ' --outReadsUnmapped Fastx '

    STAR_util.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)

    cmd = "samtools index {out_BAM}".format(out_BAM=out_BAM)
    STAR_util.logger.info(cmd)
    os.system(cmd)


def get_scope_bc(bctype):
    root_path = os.path.dirname(celescope.__file__)
    linker_f = glob.glob(f'{root_path}/data/chemistry/{bctype}/linker*')[0]
    whitelist_f = f'{root_path}/data/chemistry/{bctype}/bclist'
    return linker_f, whitelist_f


def fastq_line(name, seq, qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'


def s_common(parser):
    """subparser common arguments
    """
    parser.add_argument('--outdir', help='output dir', required=True)
    parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument('--sample', help='sample name', required=True)
    parser.add_argument('--thread', default=4)
    parser.add_argument('--debug', help='debug', action='store_true')
    return parser


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
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")


def otsu_min_support_read(array, otsu_plot):
    """
    get otsu threshold and plot
    """
    array = np.log10(array)
    hist = array2hist(array)
    thresh = threshold_otsu(hist)
    makePlot(hist, thresh, otsu_plot)
    threshold = round(10 ** thresh,1)
    return threshold
