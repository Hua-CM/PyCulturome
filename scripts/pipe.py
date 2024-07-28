# -*- coding: utf-8 -*-
# @File    :   PyCulturome.py
# @Time    :   2024/07/23 14:25:39
# @Author  :   Zhongyi Hua
# @Usage   :
# @Note    :
# @E-mail  :   njbxhzy@hotmail.com

import os
os.environ['QT_QPA_PLATFORM']='offscreen' # Set at first

import argparse
import logging
from pathlib import Path
import subprocess as sp
import sys

import pandas as pd

from assign_sample import parse_barcodes, change_seqid
from summary_library import summary_library

logger = logging.getLogger('Culturome')
logger.setLevel('INFO')
stream_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
stream_handler.setFormatter(formatter)
# 使用setHandler对记录器的处理器进行设置
logger.addHandler(stream_handler)



def parse_args():
    """Parse arguments
    """
    parser = argparse.ArgumentParser(description="Example with -1 and -2 options")
    parser.add_argument('-1', dest='fq1', required=True, type=Path,
                        help='<file_path> read1 input file path')
    parser.add_argument('-2', dest='fq2', required=True, type=Path,
                        help='<file_path> read2 input file name')
    parser.add_argument('-f', '--fwd', required=True, type=Path,
                        help='<file_path> Forward barcodes file.')
    parser.add_argument('-r', '--rev', required=True, type=Path,
                        help='<file_path> Reverse barcodes file.')
    parser.add_argument('-d', '--db', required=True, type=Path,
                        help='<file_path> The database for sintax algorithm.')
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help='<directory_path> The output directory path.')
    parser.add_argument('--bin1', type=Path, dest='usearch', default='usearch',
                        help='<file_path> The path to usearch binary (if it is not in PATH). Default: usearch')
    parser.add_argument('--bin2', '--bin', dest='vsearch', type=Path, default='vsearch',
                        help='<file_path> The path to vsearch binary (if it is not in PATH). Default: vsearch')    
    parser.add_argument('-p', '--positive', type=str, default='B12',
                        help='<str> Positive control well ID. Default: B12')
    parser.add_argument('-n', '--negative', type=str, default='A12',
                        help='<str> Negative control well ID. Default: A12')
    parser.add_argument('--threshold', type=float, default=0.95,
                        help='<float> Threshold for determining purified wells. Default: 0.95')
    args = parser.parse_args()
    return args

def run_pipeline(para_dct):
    """_summary_

    Args:
        para_dct (_type_): _description_
    """
    logger.info('Step1: Merge fastq1 and fastq2 reads')
    sp.run([para_dct['usearch'],
            '--fastq_mergepairs', para_dct['fq1'],
            '--reverse', para_dct['fq2'],
            '--fastqout', para_dct['outtmp'] / f'{para_dct["libname"]}.merged.fq'],
            check=True)

    logger.info('Step2: Assign reads to plates and wells.')
    forward_dct, reverse_dct = parse_barcodes(para_dct['fwd'], para_dct['rev'])
    cells_count_dct = change_seqid(para_dct['outtmp'] / f'{para_dct["libname"]}.merged.fq',
                                   para_dct['outtmp'] / f'{para_dct["libname"]}.changeid.fq',
                                   forward_dct,
                                   reverse_dct)
    pd.DataFrame(cells_count_dct).to_csv(para_dct['outtmp'] / f'{para_dct["libname"]}.readscount.tsv',
                                         index=False,
                                         sep='\t')

    logger.info('step3: Generate ASVs')
    sp.run([para_dct['usearch'],
            '--fastq_filter', para_dct['outtmp'] / f'{para_dct["libname"]}.changeid.fq',
            '--fastq_stripleft', '29',
            '--fastq_stripright', '24',
            '--fastq_maxee_rate', '0.01',
            '-fastq_minlen', '300', # add a length filter for now
            '--fastaout', para_dct['outtmp'] / f'{para_dct["libname"]}.filtered.fa'],
            check=True)
    sp.run([para_dct['usearch'],
           '--fastx_uniques', para_dct['outtmp'] / f'{para_dct["libname"]}.filtered.fa',
           '--fastaout', para_dct['outtmp'] / f'{para_dct["libname"]}.uniques.fa',
           '--relabel', 'Uni',
           '--minuniquesize', '10',
           '--sizeout'],
           check=True)

    sp.run([para_dct['usearch'],
           '-unoise3', para_dct['outtmp'] / f'{para_dct["libname"]}.uniques.fa',
           '-zotus', para_dct['output'] / f'{para_dct["libname"]}.rep.fa'],
           check=True)

    logger.info('step4: Generate asv count table')
    sp.run([para_dct['vsearch'],
            '-usearch_global', para_dct['outtmp'] / f'{para_dct["libname"]}.filtered.fa',
            '-db', para_dct['output'] / f'{para_dct["libname"]}.rep.fa',
            '-otutabout', para_dct['output'] / f'{para_dct["libname"]}.abun.tsv',
            #'-strand', 'both',
            '-id', '0.97',
            '-threads', '8'],
            check=True)

    logger.info('step5: Annotate ASV')
    sp.run([para_dct['vsearch'],
            '-sintax', para_dct['output'] / f'{para_dct["libname"]}.rep.fa',
            '-db', para_dct['db'],
            '-tabbedout', para_dct['outtmp'] / f'{para_dct["libname"]}.sintax',
            #'-strand', 'both',
            '-sintax_cutoff', '0.8'],
            check=True)

    logger.info('step6: Summary library results')
    summary_library(para_dct['output'] / f'{para_dct["libname"]}.abun.tsv',
                    para_dct['outtmp'] / f'{para_dct["libname"]}.sintax',
                    para_dct['output'],
                    para_dct['positive'],
                    para_dct['negative'],
                    para_dct['threshold'])


def get_common_prefix(fq_file1, fq_file2):
    """Get sample name based on fastq paired files

    Args:
        fq_file1 (str): _description_
        fq_file2 (str): _description_

    Returns:
        str: sample name
    """
    min_length = min(len(fq_file1), len(fq_file2))
    common_prefix = []

    for i in range(min_length):
        if fq_file1[i] == fq_file2[i]:
            common_prefix.append(fq_file1[i])
        else:
            break
    return ''.join(common_prefix).strip('._')


def main():
    """Main interface
    """
    args = parse_args()
    args_dict = vars(args)
    # prepare outdirectory and temporary directory
    if not args.output.exists():
        args.output.mkdir()
    args_dict['outtmp'] = args.output / 'tmp'
    if not args_dict['outtmp'].exists():
        args_dict['outtmp'].mkdir()
    args_dict['libname'] = get_common_prefix(args_dict['fq1'].stem, args_dict['fq2'].stem)
    run_pipeline(args_dict)

if __name__ == '__main__':
    main()
