# -*- coding: utf-8 -*-
# @File    :   check_sanger.py
# @Time    :   2024/09/25 14:31:54
# @Author  :   Zhongyi Hua 
# @Usage   :   compare sanger sequencing results between two sequencing results
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com
from collections import defaultdict
import re
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd


def format_alignment_custom(seq1id, seq2id, align1, align2, line_length=100):
    """Custom format function to display the alignment in the desired format with specified line length"""
    alignment_str = ""

    # Process sequences in chunks of line_length
    for i in range(0, len(align1), line_length):
        chunk1 = align1[i:i+line_length]
        chunk2 = align2[i:i+line_length]

        # Generate first sequence line
        alignment_str += f"{seq1id:20} {chunk1}\n"

        # Generate match/mismatch line
        match_line = ''.join(['|' if chunk1[j] == chunk2[j] and chunk1[j] != '-' else ' ' for j in range(len(chunk1))])
        alignment_str += f"{'':20} {match_line}\n"

        # Generate second sequence line
        alignment_str += f"{seq2id:20} {chunk2}\n"

        # Add an extra newline to separate chunks, optional
        alignment_str += "\n"

    return alignment_str


def compare_sanger(new_seq: SeqRecord, old_seq: SeqRecord):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(new_seq.seq, old_seq.seq)
    return format_alignment_custom(new_seq.id, old_seq.id, *alignments[0])


def check_sanger_main(para_dct):
    """Compare two sanger sequencing results based on
    'seqid' and 'seq'
    
    para_dict:
        old_path (Path): The old Sanger sequencing result **table** path.
        new_path (Path): The new Sanger sequencing result **table** path.
        old_rp (Bool): (reverse_complement) Does old sequence need reverse_complement?
        new_rp (Bool): (reverse_complement) Does new sequence need reverse_complement?
        out_path (Path): Output directory path
    """
    if not para_dct['out_path'].exists():
        para_dct['out_path'].mkdir()
    old_df = pd.read_table(para_dct['old_path'])
    new_df = pd.read_table(para_dct['new_path'])
    # generate old seq dct
    old_seq_dct = defaultdict(str)
    for _idx, _row in old_df.iterrows():
        old_seq_dct[_row['seqid']] = _row['seq']
    res_lst = []
    for _idx, _row in new_df.iterrows():
        new_seq = SeqRecord(id=f'{_row['seqid']}_new',
                            seq = Seq(_row['seq']).reverse_complement() if para_dct['new_rp'] else Seq(_row['seq']),
                            description='')
        if old_seq := old_seq_dct.get(_row['seqid']):
            old_seq =  SeqRecord(id=f'{_row['seqid']}_old',
                                 seq = Seq(old_seq).reverse_complement() if para_dct['old_rp'] else Seq(_row['seq']),
                                 description='')
            compare_res_text = compare_sanger(new_seq, old_seq)
            (para_dct['out_path'] / f"{_row['seqid']}.aligned.txt").write_text(
                compare_res_text
            )
            res_lst.append({'seqid': _row['seqid'],
                            'old_seq': str(old_seq.seq),
                            'new_seq': str(new_seq.seq)})
    # also output a summarized table
    tb_res = pd.DataFrame(res_lst)
    tb_res.to_csv(para_dct['out_path'] / 'compared_result.tsv', sep='\t', index=False)