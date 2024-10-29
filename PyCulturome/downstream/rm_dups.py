# -*- coding: utf-8 -*-
# @File    :   rm_dups.py
# @Time    :   2024/10/28 16:44:26
# @Author  :   Zhongyi Hua
# @Usage   :   Identify clones to isolate based on **sequence**, not OTU name
#              because some OTU have intraspeceis variance
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com

from collections import defaultdict
from itertools import combinations
import logging
from pathlib import Path
from typing import List


from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd

from .utils import runCommand

logger = logging.getLogger(__name__)

def calculate_distance(str1, str2, gap='-'):
    """_summary_

    Args:
        str1 (str): _description_
        str2 (str): _description_

    Raises:
        ValueError: _description_

    Returns:
        int : 0: Two sequences are identical; 1: Two sequences are different
    """
    if len(str1) != len(str2):
        raise ValueError("The two strings must be of the same length")

    str1 = str1.strip(gap)
    str2 = str2.strip(gap)
    if len(str1) < len(str2):
        str_short = str1
        str_long = str2
    else:
        str_short = str2
        str_long = str1
    if str_short in str_long:
        return 0
    return 1


def muscle_align(seq_lst: List[SeqRecord], tmp_dir: Path, bin_path='muscle'):
    """Align multiple sequences in a temperoary directory

    Args:
        seq_lst (List[SeqRecord]): _description_
        tmp_dir (Path): _description_
        bin_path (str, optional): _description_. Defaults to 'muscle'.

    Returns:
        Align.MultipleSeqAlignment: _description_
    """
    # using MUSCLE5 to align
    SeqIO.write(seq_lst, tmp_dir / 'input.fasta', 'fasta')
    cmd_line = [bin_path, '-align', str(tmp_dir / 'input.fasta'), '-output', str(tmp_dir / 'output.fasta')]
    runCommand(cmd_line)
    #sp.run(cmd_line, text=True, check=True)
    aligned_seqs = AlignIO.read(tmp_dir / 'output.fasta', 'fasta')
    return aligned_seqs


def rm_dup_seqs(check_rec_lst: List[SeqRecord],
                tmp_dir: Path,
                bin_path: Path,
                existing_rec_lst: List[SeqRecord] = None):
    """_summary_

    Args:
        seq_rec_lst (List[SeqRecord]): The sequences need to check
        existing_rec_lst (List[SeqRecord]): The existing SeqRecord list (if any)
        tmp_dir (Path): The temporary directory path
        bin_path (Path): The MUSCLE5 bin path

    Returns:
        seq_keep_lst : List[SeqRecord]
        seq_rm_lst   : List[SeqRecord]
    """
    if existing_rec_lst:
        seq_rec_lst = check_rec_lst + existing_rec_lst
    else:
        seq_rec_lst = check_rec_lst
    aligned_seqs = muscle_align(seq_rec_lst, tmp_dir, bin_path)
    # hard cut for now is 70 bp
    seq_dict = defaultdict(str)
    seq_rm_lst = []
    rm_id_lst = []
    for _seqrecord in aligned_seqs[:,70:-70]:
        if str(_seqrecord.seq) in seq_dict:
            rm_id_lst.append(_seqrecord.id)
        else:
            seq_dict[str(_seqrecord.seq)] = _seqrecord.id
   
    for str1, str2 in combinations(list(seq_dict.keys()), 2):
        try:
            if calculate_distance(str1, str2) == 0:
                seq_rm_lst.append(str1 if len(str1) < len(str2) else str2)
        except ValueError:
            # Handle the case where strings are of different lengths
            print(f"Skipping pair ({str1}, {str2}): Strings have different lengths")
    # filter sequences
    rm_id_lst += [seq_dict.get(_seq) for _seq in seq_rm_lst]
    seq_keep_lst = [_seq for _seq in seq_rec_lst if _seq.id not in rm_id_lst]
    seq_rm_lst = [_seq for _seq in seq_rec_lst if _seq.id in rm_id_lst]
    return seq_keep_lst, seq_rm_lst


def rm_dup_seqs2(check_rec_lst: List[SeqRecord],
                 existing_rec_lst: List[SeqRecord],
                 tmp_dir: Path,
                 bin_path: Path):
    check_id_lst = [_.id for _ in check_rec_lst]
    existing_id_lst = [_.id for _ in existing_rec_lst]
    seq_rec_lst = check_rec_lst + existing_rec_lst
    # !!IMPORTANT!!: The muscle output is out-of-order
    aligned_seqs = muscle_align(seq_rec_lst, tmp_dir, bin_path)
    seq_dict = defaultdict(str)
    check_rec_aligned_lst = []
    # construct the existing database and check data list
    for _seqrecord in aligned_seqs[:,70:-70]:
        if _seqrecord.id in existing_id_lst:
            seq_dict[str(_seqrecord.seq)] = _seqrecord.id
        else: # check
            check_rec_aligned_lst.append(_seqrecord)
    # check
    rm_id_lst = []
    for _seqrecord in check_rec_aligned_lst:
        if str(_seqrecord.seq) in seq_dict:
            rm_id_lst.append(_seqrecord.id)
        else:
            seq_dict[str(_seqrecord.seq)] = _seqrecord.id

    seq_rm_lst = []
    for str1, str2 in combinations(list(seq_dict.keys()), 2):
        try:
            if calculate_distance(str1, str2) == 0:
                # Don't remove records in existing database
                # str1 and str2 cannot exist in existing database simultaneously
                if seq_dict.get(str1) in existing_id_lst:
                    seq_rm_lst.append(str2)
                elif seq_dict.get(str2) in existing_id_lst:
                    seq_rm_lst.append(str1)
                else:
                    seq_rm_lst.append(str1 if len(str1) < len(str2) else str2)
        except ValueError:
            # Handle the case where strings are of different lengths
            print(f"Skipping pair ({str1}, {str2}): Strings have different lengths")
    # filter sequences
    rm_id_lst += [seq_dict.get(_seq) for _seq in seq_rm_lst]
    seq_keep_lst = [_seq for _seq in seq_rec_lst if _seq.id not in rm_id_lst]
    seq_rm_lst = [_seq for _seq in seq_rec_lst if _seq.id in rm_id_lst]
    # The keep_lst include existing database ids
    return seq_keep_lst, seq_rm_lst
    


def iden_isolate_main(para_dict):
    """
    NOTE: Consider iterative constructing of the bacterial culture collection
    
    para_dict:
        sanger_path (Path): The summarized top-hit Sanger sequencing result path. MUST WITH *seq* column
          There should be at least four columns in this file with column names:
            plate  : The row of cultivated 96-well plate
            well   : The column of cultivated 96-well plate
            clone  : The No. of clone for this well
            seq    : The sequence
        existing_path: An existing table with columns mentioned above and an additional 
            *ToIsolate* column without NULL.
        muscle_path (Path): MUSCLE5 bin path
        out_path (Path): Output result **directory**.
        tmp_dir (Path): temporary directory path.
    """
    tb_sanger_raw = pd.read_table(para_dict['sanger_path'])
    tb_sanger = tb_sanger_raw[['plate', 'well', 'clone', 'seq']].copy()
    # if there is seqid in row table
    tb_sanger['seqid'] = tb_sanger['plate'].astype(str) + '_' + tb_sanger['well'].astype(str) + '_' + tb_sanger['clone'].astype(str)
    # add seqid
    check_seq_lst = [SeqRecord(id=_row['seqid'], seq=Seq(_row['seq']), description='')
                     for _idx, _row in tb_sanger.iterrows()]
    if para_dict['existing_path']:
        tb_existing_raw = pd.read_table(para_dict['existing_path'])
        tb_existing_raw_used = tb_existing_raw[~(tb_existing_raw['ToIsolate'] == 0)]
        tb_existing = tb_existing_raw_used[['plate', 'well', 'clone', 'seq']].copy()
        tb_existing['seqid'] = tb_existing['plate'].astype(str) + '_' + tb_existing['well'].astype(str) + '_' + tb_existing['clone'].astype(str)
        existing_seq_lst = [SeqRecord(id=_row['seqid'], seq=Seq(_row['seq']), description='')
                            for _idx, _row in tb_existing.iterrows()]
        keep_seq_lst, rm_seq_lst = rm_dup_seqs2(check_seq_lst,
                                                existing_seq_lst,
                                                para_dict['tmp_dir'],
                                                para_dict['muscle_path'])
    else:
        # For security reasons, use 'plate_well_clone' as the identifier, disregarding the 'queryid'.
        keep_seq_lst, rm_seq_lst = rm_dup_seqs(check_seq_lst, para_dict['tmp_dir'], para_dict['muscle_path'])
    rm_seq_id_lst = [_.id for _ in rm_seq_lst]
    for _seq_id in rm_seq_id_lst:
        logger.info(f'{_seq_id} is duplicated')
    # Also output non-redundant table
    tb_sanger['ToIsolate'] = 0
    tb_sanger.loc[~tb_sanger['seqid'].isin(rm_seq_id_lst), 'ToIsolate'] = 1
    # merge to the raw table
    tb_sanger.drop(columns='seqid', inplace=True)
    tb_sanger_raw = tb_sanger_raw.merge(tb_sanger)
    # merge result and output
    if para_dict['existing_path']:
        tb_sanger_raw = pd.concat([tb_existing_raw,
                                   tb_sanger_raw])
     # Move seq to the last column
    columns = list(tb_sanger_raw.columns)
    columns.append(columns.pop(columns.index('seq')))
    tb_sanger_raw = tb_sanger_raw[columns]
    tb_sanger_raw.to_csv(
        para_dict['out_path'] / 'Isolation.tsv',
        sep='\t',
        index=False
    )
    SeqIO.write(
        keep_seq_lst,
        para_dict['out_path'] / 'Bacteria_ToIsolate.fasta',
        'fasta')

    