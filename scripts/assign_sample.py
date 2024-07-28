# -*- coding: utf-8 -*-
# @Time    : 2024/7/12
# @Author  : Zhongyi Hua
# @File    : assign_sample.py
# @Usage   : assign plates and cells name to the sequence id
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com
import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict


from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


def deter_ambiguous_barcode(barcode, barcode_dct):
    """Determine ambiguous barcode using reverse match

    Args:
        barcode (_type_): _description_
        barcode_dct (_type_): _description_
        direct (_type_): _description_
    """
    _idx = barcode.find('N')
    tmp_barcode = barcode[:_idx] + barcode[_idx+1:]
    tmp_barcode_dct = {_key[:_idx]+_key[_idx+1:]: _val for _key, _val in barcode_dct.items()}
    tmp_id = tmp_barcode_dct.get(tmp_barcode)
    if not tmp_id:
        tmp_id = 'Unclassified'
    return tmp_id


def change_seqid(merged_fq_path: Path, output_fq_path: Path, fwd_bar_dct: Dict, rev_bar_dct: Dict):
    """Change sequence id in FASTQ file to 'sample_id.<number>' format for subsequent analysis.

    Args:
        merged_fq_path (Path): The merged FASTQ file path (by usearch) 
        output_fq_path (Path): The output FASTQ file path (for search)
        fwd_bar_dct (Dict): Generated from parse_barcodes
        rev_bar_dct (Dict): Generated from parse_barcodes

    Returns:
        sample_counter: {'<cell number>_<plate number>': <reads_counter>, ...}. 
                        e.g. {'F1_P1': 12345}: 
    """
    merged_amps=SeqIO.parse(merged_fq_path, 'fastq')
    sample_counter = defaultdict(lambda: 1)
    out_seq_lst = []

    fwd_bar_len = len(list(fwd_bar_dct.keys())[0])
    rev_bar_len = len(list(rev_bar_dct.keys())[0])

    for _seq in merged_amps:
        cellid = fwd_bar_dct.get(str(_seq.seq[0: fwd_bar_len]))
        plateid = rev_bar_dct.get(str(_seq.seq[-rev_bar_len:]))
        if not cellid:
            cellid = deter_ambiguous_barcode(str(_seq.seq[0: fwd_bar_len]), fwd_bar_dct)
        if not plateid:
            plateid = deter_ambiguous_barcode(str(_seq.seq[-rev_bar_len:]), rev_bar_dct)
        sampleid = f'{cellid}_{plateid}'
        sample_seq_num = sample_counter.get(sampleid, 1)
        _seq.id = f'{sampleid}.{sample_seq_num}'
        _seq.description = ''
        out_seq_lst.append(_seq)
        sample_counter[sampleid] += 1
    SeqIO.write(out_seq_lst, output_fq_path, 'fastq')
    sample_count_lst = []
    for _key, _val in sample_counter.items():
        _cell, _plate = _key.split('_')
        sample_count_lst.append({'cell':_cell, 'plate': _plate, 'reads': _val})
    return sample_count_lst


def parse_barcodes(fwd_bar_path: Path, rev_bar_path: Path):
    """Parse barcodes to plates and cells

    Args:
        fwd_bar_path (Path): Forward barcode path
        rev_bar_path (Path): Reverse barcode path
    Returns:
        bar_id_dct (Dict): 
    """
    def read_barcode_file(file_path):
        return {_seq: _id for _id, _seq in
                (line.split('\t') for line in file_path.read_text().strip().split('\n') if line)
               }

    fwd_dct = read_barcode_file(fwd_bar_path)
    rev_dct = read_barcode_file(rev_bar_path)
    # change seqs in rev_dct to reverse complement
    rev_dct = {str(Seq(_seq).reverse_complement()): _id
               for _seq, _id in rev_dct.items()}

    #bar_id_dct = {(k1, k2): f'{v1}_{v2}' for (k1, v1), (k2, v2)
    #  in product(fwd_dct.items(), rev_dct.items())}
    return fwd_dct, rev_dct


def parse_args():
    """Parse CLI commands
    """
    parser = argparse.ArgumentParser(description="Assign plates and cells name to the sequence \
                                     id for subsequence analysis")
    parser.add_argument('-i', '--input', required=True, type=Path,
                        help='<file_path>  The input merged FASTQ file.')
    parser.add_argument('-f', '--fwd', required=True, type=Path,
                        help='<file_path> Forward barcodes file.')
    parser.add_argument('-r', '--rev', required=True, type=Path,
                        help='<file_path> Reverse barcodes file.')
    parser.add_argument('-o', '--output1', required=True, type=Path,
                        help='<file_path> The output FASTQ file.')
    parser.add_argument('-O', '--output2', required=True, type=Path,
                        help='<file_path> The output reads count.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args_global = parse_args()
    forward_dct, reverse_dct = parse_barcodes(args_global.fwd, args_global.rev)
    cells_count_dct = change_seqid(args_global.input,
                                   args_global.output1,
                                   forward_dct,
                                   reverse_dct)
    pd.DataFrame(cells_count_dct).to_csv(args_global.output2, index=False, sep='\t')
