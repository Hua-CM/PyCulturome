# -*- coding: utf-8 -*-
# @File    :   rm_dups.py
# @Time    :   2024/09/09 18:11:06
# @Author  :   Zhongyi Hua
# @Usage   :   
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com

from collections import defaultdict
from itertools import combinations
import logging
from pathlib import Path
import subprocess as sp
import sys
from typing import List


from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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


def runCommand(cmd, timeout=None):
    """ run shell command

    @param cmd: command to execute
    @param timeout: timeout for command execution

    @return: (return code from command, command output)
    """
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)
    #output = ''
    for line in p.stdout:
        line = line.decode(encoding='utf-8').strip()
        for _line in line.split('\n'):
            sys.stdout.write(_line + '\n')
            sys.stdout.flush()
    p.wait(timeout)


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


def rm_dup_seqs(seq_rec_lst, tmp_dir, bin_path):
    """_summary_

    Args:
        seq_lst (List[SeqRecord]): _description_

    Returns:
        _type_: _description_
    """
    aligned_seqs = muscle_align(seq_rec_lst, tmp_dir, bin_path)
    # hard cut for now 70 bp
    seq_dict = defaultdict(str)
    rm_seq_lst = []
    rm_id_lst = []
    for _seqrecord in aligned_seqs[:,70:-70]:
        if str(_seqrecord.seq) in seq_dict:
            rm_id_lst.append(_seqrecord.id)
        else:
            seq_dict[str(_seqrecord.seq)] = _seqrecord.id
   
    for str1, str2 in combinations(list(seq_dict.keys()), 2):
        try:
            if calculate_distance(str1, str2) == 0:
                rm_seq_lst.append(str1 if len(str1) < len(str2) else str2)
        except ValueError:
            # Handle the case where strings are of different lengths
            print(f"Skipping pair ({str1}, {str2}): Strings have different lengths")
    # filter sequences
    rm_id_lst += [seq_dict.get(_seq) for _seq in rm_seq_lst]
    seq_keep_lst = [_seq for _seq in seq_rec_lst if _seq.id not in rm_id_lst]
    return seq_keep_lst
