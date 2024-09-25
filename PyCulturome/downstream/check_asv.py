# -*- coding: utf-8 -*-
# @File    :   check_asv.py
# @Time    :   2024/09/23 17:17:48
# @Author  :   Zhongyi Hua 
# @Usage   :   Evaluate whether sequences that do not align with the reference 
#              ASV can be reassigned to a different ASV.
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com

import logging
from pathlib import Path

from Bio import SeqIO
import pandas as pd

from .utils import runCommand

logger = logging.getLogger(__name__)


def run_blastn(query_path: Path, db_path: Path, out_path: Path, para_dict: dict):
    blast_cmd = [para_dict['bin'],
                 '-query', str(query_path),
                 '-db', str(db_path),
                 '-out', str(out_path)]
    for _key, _val in para_dict.items():
        if _key == 'bin':
            continue
        blast_cmd.append(f'-{_key}')
        blast_cmd.append(str(_val))
    runCommand(blast_cmd)


def run_makeblastdb(in_path: Path, out_path: Path, bin_path):
    blast_cmd = [str(bin_path), '-in', str(in_path), '-dbtype', 'nucl', '-out', str(out_path)]
    runCommand(blast_cmd)


def check_asv_main(para_dict):
    """_summary_

    para_dict:
        query_path (Path): Path to query sequences fasta
        asv_path (Path): Path to ASV fasta
        db_path (Path): Path to NCBI 16S rRNA gene database
        out_path (Path): Result table path
        blast_dir: BLAST software directory
        tmp_dir: temporary directory
        threads: Number of threads
    """
    if not para_dict['tmp_dir'].exists():
        para_dict['tmp_dir'].mkdir()

    run_makeblastdb(para_dict['asv_path'],
                    para_dict['tmp_dir'] / 'asv_db',
                    str(para_dict['blast_dir'] / 'makeblastdb'))
    
    blastn_para_dct = {
        'bin': str(para_dict['blast_dir'] / 'blastn'),
        'evalue': 1e-5,
        'outfmt': '6 qacc sacc slen length pident evalue',
        'max_hsps': 1,
        'max_target_seqs': 1,
        'num_threads': str(para_dict['threads']),
    }
    
    run_blastn(para_dict['query_path'],
               para_dict['tmp_dir'] / 'asv_db',
               para_dict['tmp_dir'] / 'asv_res.tsv',
               blastn_para_dct)
    # parse query ASV result
    tb_res_asv= pd.read_table(para_dict['tmp_dir'] / 'asv_res.tsv',
                              names=['queryid', 'subjectid', 'slen', 'alignlen',
                                     'identity', 'evalue'])
    res_lst = []
    query_ncbi_lst = []
    for _idx, _row in tb_res_asv.iterrows():
        if (_row['slen'] == _row['alignlen']) and (_row['identity'] == 100):
            res_lst.append({'seqid': _row['queryid'],
                            'subjectid': _row['subjectid']})
        else:
            query_ncbi_lst.append(_row['queryid'])
    # Query NCBI
    query_seqs = SeqIO.to_dict(SeqIO.parse(para_dict['query_path'], 'fasta'))
    query_seqs = [query_seqs.get(_) for _ in query_ncbi_lst]
    SeqIO.write(query_seqs, para_dict['tmp_dir'] / 'query.fasta', 'fasta')
    blastn_para_dct2 = {
        'bin': str(para_dict['blast_dir'] / 'blastn'),
        'evalue': 1e-5,
        'outfmt': "6 qacc sacc qlen slen length qstart qend sstart send ssciname staxid pident evalue",
        'max_hsps': 1,
        'max_target_seqs': 1,
        'num_threads': str(para_dict['threads']),
    }    
    run_blastn(para_dict['tmp_dir'] / 'query.fasta',
               para_dict['db_path'],
               para_dict['tmp_dir'] / 'ncbi_res.tsv',
               blastn_para_dct2)
    tb_res_ncbi= pd.read_table(para_dict['tmp_dir'] / 'ncbi_res.tsv',
                               names=['queryid', 'subjectid', 'qlen',   'slen', 'alignlen',
                                      'qstart',  'qend',      'sstart', 'send', 'species',
                                      'taxid', 'identity', 'evalue'])
    for _idx, _row in tb_res_ncbi.iterrows():
        res_lst.append({'seqid': _row['queryid'],
                        'subjectid': _row['subjectid'],
                        'species': _row['species'],
                        'taxid': _row['taxid'],
                        'identity': _row['identity']})
    # merge ASV and NCBI result and output
    tb_res = pd.DataFrame(res_lst)
    tb_res.to_csv(para_dict['out_path'], sep='\t', index=False)