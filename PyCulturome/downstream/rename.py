# -*- coding: utf-8 -*-
# @File    :   rename.py
# @Time    :   2024/09/10 16:32:28
# @Author  :   Zhongyi Hua
# @Usage   :   Renaming Strains for Long-Term Preservation
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com
import logging

from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm

import pandas as pd

logger = logging.getLogger(__name__)

def change_seq_id(seq_path: Path, meta_path: Path):
    seq_dct = SeqIO.to_dict(SeqIO.parse(seq_path, 'fasta'))

    tb_info = pd.read_table(meta_path)
    info_lst = tb_info.to_dict('records')

    out_lst = []

    for _info_dct in tqdm(info_lst):
        _newid = _info_dct.get('newid')
        _seq = seq_dct.get(_info_dct.get('oldid'))
        _seq.id = _newid
        _seq.description = ''
        out_lst.append(_seq)
    return out_lst


def rename_main(para_dct):
    """_summary_

    para_dct:
        in_path (Path): The raw sequence fasta path.
        meta_path (Path): A TSV file contanins 'oldid' and 'newid'.
        out_path (Path): The converted sequence fasta path.
    """
    logger.info('Start changing sequences id')
    out_seq_lst = change_seq_id(para_dct['in_path'], para_dct['meta_path'])
    SeqIO.write(out_seq_lst, para_dct['out_path'], 'fasta')
    logger.info('Sequences id have been changed!')
