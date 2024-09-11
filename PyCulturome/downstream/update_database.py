# -*- coding: utf-8 -*-
# @File    :   update_database.py
# @Time    :   2024/09/10 09:43:43
# @Author  :   Zhongyi Hua 
# @Usage   :   
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com
import argparse
from collections import defaultdict
from itertools import product
import logging
from pathlib import Path
import sys

from Bio import SeqIO
import pandas as pd

from .rm_dups import calculate_distance,  muscle_align

logger = logging.getLogger(__name__)

def updata_database(db_path, new_seq_path, tmp_dir, bin_path):
    """_summary_

    Args:
        db_path (_type_): _description_
        new_seq_path (_type_): _description_
        tmp_dir (_type_): _description_
        bin_path (_type_): _description_

    Returns:
        List[SeqRecord]: _description_
    """
    seqrec_existing = list(SeqIO.parse(db_path, 'fasta'))
    seqrec_adding = list(SeqIO.parse(new_seq_path, 'fasta'))
    exist_id_lst = [_.id for _ in seqrec_existing]
    aligned_seqs = muscle_align((seqrec_existing + seqrec_adding),
                                tmp_dir,
                                bin_path)
    seq_existing = defaultdict(str)
    seq_adding = defaultdict(str)
    for _seq in aligned_seqs[:,70:-70]:
        if _seq.id in exist_id_lst:
            seq_existing[str(_seq.seq)] = _seq.id
        else:
            seq_adding[str(_seq.seq)] = _seq.id
    
    # remove duplicated sequences in adding sequences first
    seq_dct = defaultdict(str)
    rm_id_lst = []
    for _seq, _seqid in seq_adding.items():
        if _seq[70:-70] in seq_dct:
            rm_id_lst.append(_seqid)
        else:
            seq_dct[_seq] = _seqid
    
    for str1, str2 in product(list(seq_existing.keys()), list(seq_dct.keys())):
        if calculate_distance(str1, str2) == 0:
            rm_id_lst.append(seq_dct.get(str2))
    
    seq_add_lst = [_seq for _seq in seqrec_adding if _seq.id not in rm_id_lst]
    return seq_add_lst


def parse_args():
    """Parse arguments
    """
    parser = argparse.ArgumentParser(description="Example with -1 and -2 options")
    parser.add_argument('-i', '--input', required=True, type=Path, dest='in_path',
                        help='<file_path> The sequences to be added')
    parser.add_argument('-d', '--db', required=True, type=Path, dest='db_path',
                        help='<file_path> The existing 16S rRNA gene database. fasta file')
    parser.add_argument('-o', '--output', required=True, type=Path, dest='out_path',
                        help='<file_path> The output fasta path.')
    parser.add_argument('-m', '--meta', type=Path, default=None, dest='meta_path',
                        help='<file_path> The meta information for sequences to be added (if any).'
                        'Note: there must be a column named "seqid" conating sequences id')
    parser.add_argument('--bin', type=Path, dest='bin_path', default=Path(''),
                        help='<file_path> The path to MUSCLE5 binary (if it is not "muscle" in PATH)') 
    parser.add_argument('--tmp', type=Path, dest='tmp_dir', default=Path('tmp'),
                        help='<int> Temperory directory path. Default: ./tmp')
    args = parser.parse_args()
    return args


def update_db_main(para_dict):
    """_summary_

    Args:
        db_path (Path): _description_
        in_path (Path): _description_
        out_path (Path): _description_
        tmp_dir (Path): _description_
        bin_path (Path): _description_
        meta_path (Path): _description_
    """
    if not para_dict['tmp_dir'].exists():
        para_dict['tmp_dir'].mkdir()

    logger.info('Database update start, please wait for a moment')
    seq_add_lst = updata_database(para_dict['db_path'],
                                  para_dict['in_path'],
                                  para_dict['tmp_dir'],
                                  para_dict['bin_path'])
    SeqIO.write(seq_add_lst, para_dict['out_path'], 'fasta')
    seq_id_lst = [_seq.id for _seq in seq_add_lst]
    if para_dict['meta_path']:
        logger.info('Filter meta information')
        tb_meta = pd.read_table(para_dict['meta_path'])
        tb_meta = tb_meta[tb_meta['seqid'].isin(seq_id_lst)]
        out_path = para_dict['out_path'].parent / (para_dict['out_path'].stem + '.tsv')
        tb_meta.to_csv(out_path, sep='\t', index=False)
        logger.info('Filter meta information done')


if __name__ == '__main__':
    ARGS = parse_args()
    update_db_main(vars(ARGS))
