# -*- coding: utf-8 -*-
# @File    :   mergesanger2.py
# @Time    :   2024/09/03 09:57:51
# @Author  :   Zhongyi Hua
# @Usage   :   Merge two sanger sequencing reads to the contig
# @Note    :   Modified from https://github.com/shiqiang-lin/merge_sanger_sequences_v3
# @E-mail  :   njbxhzy@hotmail.com

from pathlib import Path
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_seq(file_path: Path):
    """Read sequence from file. All of fasta/seq/ab1 is OK.

    Args:
        file_path (Path): The read file path

    Returns:
        str: The character of sequence
    """
    match file_path.suffix:
        case '.ab1':
            seq_str = str(SeqIO.read(file_path, 'abi').seq)
        case '.fasta' | 'fa':
            seq_str = str(SeqIO.read(file_path, 'fasta').seq)
        case '.seq':
            seq_str = file_path.read_text().strip().replace('\n', '')
    return seq_str


def pairwise_align(seq1: str, seq2: str, gap_open=10, gap_extend=5):
    """
    align two sequences using Align.PairwiseAligner() function.
        Be mindful of these parameter settings
    
    Args:
        seq1 (str): _description_
        seq2 (str): _description_
        gap_open (_type_): _description_
        gap_extend (_type_): _description_

    Returns:
        Bio.Align.Alignment: _description_
    """
    # Create a PairwiseAligner object
    aligner = Align.PairwiseAligner()

    # Set parameters
    aligner.mode = 'global'
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -gap_open
    aligner.extend_gap_score = -gap_extend
    aligner.end_gap_score = 0
    # Perform alignment
    alignments = aligner.align(seq1, seq2)

    best_alignment = alignments[0]
    return best_alignment


def find_n_consecutive(aligned_seq1, aligned_seq2, n):
    """
    find the position of n consecutive matching fragment
    """
    if len(aligned_seq1) != len(aligned_seq2):
        print("length of aligned seq1 must equal to aligned seq2.")
        return -1
    for i in range(len(aligned_seq1) - n + 1):
        if aligned_seq1[i:i + n] == aligned_seq2[i:i + n]:
            return i
    return -1


def merge_two_seq(seq1, seq2):
    """align seq1 and seq2
    Args:
        seq1 (str): _description_
        seq2 (str): _description_

    Returns:
        str: _description_
    """
    # Perform Needleman-Wunsch alignment
    alignment = pairwise_align(seq1, seq2)
    # Format alignment for output

    #get aligned seq1 and seq2
    aligned_fasta = alignment.format("fasta").split('\n')
    aligned_seq1 = aligned_fasta[1].strip()
    aligned_seq2 = aligned_fasta[3].strip()

    consensus_len = 50 # the consensus base number set to 50 for now
    aligned_pos = find_n_consecutive(aligned_seq1, aligned_seq2, consensus_len)
    
    if aligned_pos == -1:
        print("can not find consecutive %d equal in two sequences.", consensus_len)
    
    # Determine order
    seq1_num = aligned_seq1[:aligned_pos].count('-')
    seq2_num = aligned_seq2[:aligned_pos].count('-')
    if seq1_num > seq2_num: # seq2 in front
        merged_seq_str = aligned_seq2[:aligned_pos+consensus_len].replace('-', '') + aligned_seq1[aligned_pos+consensus_len:].replace('-', '')
    else:  # seq1 in front
        merged_seq_str = aligned_seq1[:aligned_pos+consensus_len].replace('-', '') + aligned_seq2[aligned_pos+consensus_len:].replace('-', '')    
    return merged_seq_str


def merge_sanger(seq1_path: Path, seq2_path: Path, seq_id='merged'):
    """merge Sanger sequencing result.
    
    **IMPORTANCE**: This function always treats the 
    direction of read1 as the forward orientation.

    Args:
        seq1_path (Path): The read1 path
        seq2_path (Path): The read2 path
    """
    seq1_str = read_seq(seq1_path)
    seq2_str = read_seq(seq2_path)
    seq2r_str = str(Seq(seq2_str).reverse_complement())
    merged_seq_str = merge_two_seq(seq1_str, seq2r_str)
    merged_seq = SeqRecord(seq=Seq(merged_seq_str),
                           id=seq_id,
                           name='',
                           description='')
    return merged_seq
