import argparse
import json
import os
from collections import defaultdict

import pandas as pd
from celescope.tools import utils


@utils.add_log
def append_sample_data(sample, sample_data_dict, all_data_dict, steps):
    """
    append sample data to all_data_dict
    Args:
        all_data_dict - key: step, value: list of sample series
    """
    for step in steps:
        title = step + '_summary'
        if title in sample_data_dict:
            metric_list = sample_data_dict[title]['metric_list']
            if metric_list:
                metric_dict = {'sample': sample}
                for metric in metric_list:
                    name = metric['name'].replace(' ', '_')
                    metric_dict[name] = metric['display']

                all_data_dict[step].append(pd.Series(metric_dict))

@utils.add_log
def write_merge_report(all_data_dict, merge_report_handle, steps):
    """
    merge_report:
        ## sample_summary
        sample,sample_ID,Assay...
        sample_1,sample_1,rna...
        ## barcode_summary

    Args:
        all_data_dict - key: title, value: list of sample series
    """

    for step in steps:
        if all_data_dict[step]:
            merge_report_handle.write(f'## {step}_summary\n')
            df = pd.concat(all_data_dict[step], axis=1).T
            df.to_csv(merge_report_handle, index=False, mode='a', sep='\t')
            merge_report_handle.write('\n')

@utils.add_log
def run(args):
    # out
    os.chdir(args.outdir)
    out_file = 'merge.xls'

    all_data_dict = defaultdict(list) 
    steps = args.steps.split(',')
    samples = args.samples.split(',')
    run.logger.info(f'samples: {samples}')

    for sample in samples:
        sample_data_dict = json.load(open(f'{sample}/.data.json', 'r'))
        append_sample_data(sample, sample_data_dict, all_data_dict, steps)

    write_merge_report(all_data_dict, open(out_file, 'w'), steps)


def main():
    parser = argparse.ArgumentParser('merge report')
    parser.add_argument('--outdir', help='outdir', required=True)
    parser.add_argument('--samples', help='samples, seperated by comma', required=True)
    parser.add_argument('--steps', help='steps', required=True)
    parser.add_argument('--rm_files', action='store_true', help='remove all fq and bam after running')
    args = parser.parse_args()

    run(args)

@utils.add_log
def rm_files():
    cmd = '''
        find . -iname '*.fq*' -delete;
        find . -iname '*.bam' -not -path './*/*.featureCounts/*name_sorted.bam' -delete;
    '''
    os.system(cmd)


if __name__ == '__main__':
    main()
