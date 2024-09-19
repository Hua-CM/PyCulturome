# -*- coding: utf-8 -*-
# @File    :   validate_asv.py
# @Time    :   2024/09/18 21:23:27
# @Author  :   Zhongyi Hua 
# @Usage   :   Validate whether the isolated strain is identical to the sequenced ASV
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com
import logging
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio import SeqIO
import pandas as pd

logger = logging.getLogger(__name__)

def pairwise_align(seq1, seq2):
    aligner = PairwiseAligner()
    # 参数设置
    aligner.mode = 'global'  # 16S rRNA基因比对通常使用全局比对
    aligner.open_gap_score = -5  # 较高的gap打开罚分，这样更严格
    aligner.extend_gap_score = -2  # 较高的gap延伸罚分
    aligner.match_score = 2  # 较高的匹配得分
    aligner.mismatch_score = -2  # 较高的不匹配罚分
    alignments = aligner.align(seq1, seq2)
    return alignments[0]


def check_identical(aligned_seq1, aligned_seq2, identity: float = 0.99, gap='-'):
    """check whether the isolated strain sequence identical to the 

    Args:
        aligned_seq1 (str): The aligned sequence1
        aligned_seq2 (str): The aligned sequence2
        identity (float, optional): _description_. Defaults to 0.99.
    
    Return:
        0/1: 0 represents consistency, 1 represents inconsistency
    """
    if len(aligned_seq1) != len(aligned_seq2):
        raise ValueError("The two strings must be of the same length")

    distance = 0
    for char1, char2 in zip(aligned_seq1, aligned_seq2):
        if char1 != char2 and char1 != gap and char2 != gap:
            distance += 1
    if (1- distance / len(aligned_seq1)) < identity:
        return 0
    return 1


def validate_asv(para_dict):
    """Validate whether the 
    
    para_dict:
        sanger_path (Path): The Sanger sequencing path.
        ngs_path (Path): The NGS analysis result: purified_well.tsv.
        seq_path (Path): The *.rep.fa.
        out_path (Path): The result output path.
    """
    tb_sanger = pd.read_table(para_dict['sanger_path'])
    tb_ngs = pd.read_table(para_dict['ngs_path'])
    asv_seq_dct = SeqIO.to_dict(SeqIO.parse(para_dict['seq_path'], 'fasta'))
    well_dct = {}
    for _idx, _row in tb_ngs.iterrows():
        well_dct.setdefault(_row['plate'] + _row['well'], _row['OTUID'])
    check_res_lst = []
    for _idx, _row in tb_sanger.iterrows():
        raw_well = _row['plate'] + _row['well']
        try:
            seq_sanger = _row['seq']
            seq_ngs = asv_seq_dct.get(well_dct.get(raw_well))
            alignment = pairwise_align(seq_sanger, seq_ngs)
            check_res_lst.append({'plate': _row['plate'],
                                  'well' : _row['well'],
                                  'consistency': check_identical(alignment[0], alignment[1])})
        except:
            logger.info(f'{raw_well} not in NGS result, please check')
    tb_validate = pd.DataFrame(check_res_lst)
    tb_sanger = tb_sanger.merge(tb_validate)
    # move seq to the last
    columns = list(tb_sanger.columns)
    columns.append(columns.pop(columns.index('seq')))
    tb_sanger = tb_sanger[columns]
    tb_sanger.to_csv(para_dict['out_path'], sep='\t', index=False)