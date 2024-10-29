# -*- coding: utf-8 -*-
# @File    :   sanger.py
# @Time    :   2024/09/02 22:35:01
# @Author  :   Zhongyi Hua
# @Usage   :   Parse Sanger sequencing result for 96-well plates
# @Note    :
# @E-mail  :   njbxhzy@hotmail.com
from collections import defaultdict
import logging
from pathlib import Path
import re

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
from tqdm import tqdm

from .merge_sanger import merge_sanger, read_seq
from .rm_dups import rm_dup_seqs

logger = logging.getLogger(__name__)
pd.options.display.float_format = '{:,.0f}'.format

def merge_sanger_reads(input_dir: Path, extra_str, f_primer:str, r_primer:str):
    """Merge all sanger sequencing result files in a directory

    Args:
        input_dir (Path): The folder containing sequencing files.
        extra_str (str) : The extra str need to be removed.
        f_primer  (str) : Forward primer name.
        r_primer  (str) : Reverse primer name.

    Returns:
        List[SeqRecord]: The merged seq record
    """
    sample_dct = defaultdict(defaultdict)
    for _file in input_dir.glob('*.ab1'):
        _sample = re.match(r'.*?\(', _file.stem).group(0)[:-2].replace(extra_str, '')
        if re.search(f_primer, _file.stem):
            sample_dct[_sample]['F'] = _file
        if re.search(r_primer, _file.stem):
            sample_dct[_sample]['R'] = _file
    # Merge
    seq_lst = []
    for _sample, _file_dct in tqdm(sample_dct.items()):
        _forward_file = _file_dct.get('F')
        _reverse_file = _file_dct.get('R')
        if _forward_file and _reverse_file:
            seq_lst.append(merge_sanger(_forward_file, _reverse_file, seq_id=_sample))
        elif _forward_file:
            logger.info('The reverse read of sample %s is missing, skip', _sample)
        else:
            logger.info('The forward read of sample %s is missing, skip', _sample)
    return seq_lst


def summary_sanger_reads(input_dir: Path, extra_str: str):
    seq_lst = []
    for _file in input_dir.glob('*.ab1'):
        _sample = re.match(r'.*?\(', _file.stem).group(0)[:-2].replace(extra_str, '')
        seq_str = read_seq(_file)
        seq_lst.append(SeqRecord(id = _sample, seq=Seq(seq_str), description='', name=''))
    return seq_lst


def sanger_main(para_dict):
    """_summary_

    para_dict:
        input_dir (Path): Directory containing Sanger sequencing reads end with '.ab1'
        is_merge (Bool): Should Sanger reads be merged?
        extra_str (str): Extra characters in file names that need to be removed.
        meta_path (Path): The meta info to match PCR plate and picked clones.
          There should be at least four columns in this file with column names:
            sample: The sample id
            row   : The row of PCR plate (If there not 'sample')
            column: The column of PCR plate (If there not 'sample')
            plate : The No. of plate for bacteria isolation
            well  : The No. of well for bacteria isolation
            clone : The No. of clone for this well
        keep_raw(Bool): if use 'sample' as ID, then set this para as True, otherwise False.
        db_path (Path): NCBI 16S rRNA gene database.
        out_path (Path): Output result **directory**.
        f_primer (str): Forward primer name. Default: 27F
        r_primer (str): Reverse primer name. Default: 1492R
        blast_dir (Path): The BLAST bin path
        tmp_dir (Path): temporary directory path.
        threads (int): Threads for BLAST.
    """
    if not para_dict['tmp_dir'].exists():
        para_dict['tmp_dir'].mkdir()
    
    if not para_dict['out_path'].exists():
        para_dict['out_path'].mkdir()
    if para_dict['is_merge'] == True:
        sample_seq_lst = merge_sanger_reads(para_dict['input_dir'], para_dict['extra_str'], para_dict['f_primer'], para_dict['r_primer'])
    else:
        sample_seq_lst = summary_sanger_reads(para_dict['input_dir'], para_dict['extra_str'])
    query_seq_path = Path(para_dict['tmp_dir'] / 'all_sample_merged.fa')
    SeqIO.write(sample_seq_lst, query_seq_path,'fasta')
    tb_seq = pd.DataFrame([{'queryid': _seq.id, 'seq': _seq.seq} for _seq in sample_seq_lst])
    tb_seq['queryid']= tb_seq['queryid'].astype(str)
    logger.info('Merge Sanger sequencing reads. Done!')
    # Blast
    blast_cmd = NcbiblastnCommandline(
        cmd = para_dict['blast_dir'] / 'blastn',
        query = query_seq_path,
        db=para_dict['db_path'],
        evalue=1e-5,
        outfmt="6 qacc sacc qlen slen length ssciname staxid pident evalue",
        max_hsps=1,
        max_target_seqs=20,
        num_threads=para_dict['threads'],
        out=para_dict['tmp_dir'] / 'blast.out'
    )
    blast_cmd()
    tb_blast= pd.read_table(para_dict['tmp_dir'] / 'blast.out',
                  names=['queryid', 'subjectid', 'qlen', 'slen', 'alignlen',
                         'species', 'taxid', 'identity', 'evalue'])
    tb_blast['queryid']= tb_blast['queryid'].astype(str)
    logger.info('BLAST against NCBI 16S rRNA gene database. Done!')
    # parse meta information
    tb_meta = pd.read_table(para_dict['meta_path'])
    if not para_dict['keep_raw']:
        tb_meta['queryid'] = tb_meta['row'] + tb_meta['column'].astype(str)
    else:
        tb_meta['queryid'] = tb_meta['sample'].astype(str)
    # merge
    tb_out = tb_seq.merge(tb_meta, how='left').merge(tb_blast, how='left')
    columns = list(tb_out.columns)
    columns.append(columns.pop(columns.index('seq')))
    tb_out = tb_out[columns]
    # Output
    # #Output all result
    tb_out.to_csv(para_dict['out_path'] / 'all_hit.tsv', sep='\t', index=False)
    logger.info('Parse BLAST results. Done!')
    # # Output non-redundent top results
    # ## Select the top hit for each sequence
    tb_dedup = tb_out.sort_values(by='evalue').drop_duplicates(subset=['queryid'])
    tb_dedup['seqid'] = tb_dedup['plate'].astype(str) + '_' + tb_dedup['well'].astype(str) + '_' + tb_dedup['clone'].astype(str)
    # construct the list of sequences
    seq_lst = [SeqRecord(id=_row['seqid'],
                    seq=Seq(_row['seq']),
                    description=''
                    ) for _idx, _row in tb_dedup.iterrows()]
    # Move the function of removing duplicated isolates to an independent module. 
    """
    keep_seq_lst, rm_seq_lst = rm_dup_seqs(seq_lst, para_dict['tmp_dir'], para_dict['muscle_path'])
    rm_seq_id_lst = [_.id for _ in rm_seq_lst]
    for _seq_id in rm_seq_id_lst:
        logger.info(f'{_seq_id} is duplicated')
    # Also output non-redundant table
    tb_dedup['ToIsolate'] = 0
    tb_dedup.loc[~tb_dedup['seqid'].isin(rm_seq_id_lst), 'ToIsolate'] = 1
    """
    columns = list(tb_dedup.columns)
    columns.append(columns.pop(columns.index('seq'))) # Move seq to the last column
    tb_dedup = tb_dedup[columns]
    tb_dedup.to_csv(para_dict['out_path'] / 'top_hit.tsv',
                    sep='\t',
                    index=False)
    logger.info('Summary results. Done!')
    SeqIO.write(
        seq_lst,
        para_dict['out_path'] / 'All_sequences.fasta',
        'fasta'
    )
    """
    SeqIO.write(
        keep_seq_lst,
        para_dict['out_path'] / 'Bacteria_ToIsolate.fasta',
        'fasta')
    """
    logger.info('Output fasta sequence. Done!')
    #shutil.rmtree(para_dict['tmp_dir'])
