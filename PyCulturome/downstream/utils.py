# -*- coding: utf-8 -*-
# @File    :   utils.py
# @Time    :   2024/09/23 17:26:47
# @Author  :   Zhongyi Hua 
# @Usage   :   
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com
import subprocess as sp
import sys

from Bio.Align import PairwiseAligner


def runCommand(cmd, timeout=None):
    """ run shell command

    @param cmd: command to execute
    @param timeout: timeout for command execution

    @return: (return code from command, command output)
    """
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)
    #output = ''
    for line in p.stdout:
        try:
            line = line.decode(encoding='utf-8').strip()
        except UnicodeDecodeError:
            line = line.decode(encoding='iso-8859-1').strip()
        for _line in line.split('\n'):
            sys.stdout.write(_line + '\n')
            sys.stdout.flush()
    p.wait(timeout)

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