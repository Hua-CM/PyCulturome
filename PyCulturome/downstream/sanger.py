# -*- coding: utf-8 -*-
# @File    :   sanger.py
# @Time    :   2024/09/02 22:35:01
# @Author  :   Zhongyi Hua
# @Usage   :   Parse Sanger sequencing result for 96-well plates
# @Note    :
# @E-mail  :   njbxhzy@hotmail.com
import argparse
from collections import defaultdict
import logging
from pathlib import Path
import re
import sys

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
from tqdm import tqdm

from .merge_sanger import merge_sanger
from .rm_dups import rm_dup_seqs

logger = logging.getLogger(__name__)

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


def sanger_main(para_dict):
    """_summary_

    para_dict:
        input_dir (Path): Directory containing Sanger sequencing reads end with '.ab1'
        extra_str (str): Extra characters in file names that need to be removed.
        meta_path (Path): The meta info to match PCR plate and picked clones.
          There should be five columns in this file **without column names:**
          row   : The row of PCR plate
          column: The column of PCR plate
          plate : The No. of plate for bacteria isolation
          well  : The No. of well for bacteria isolation
          clone : The No. of clone for this well
        db_path (Path): _description_
        out_path (Path): _description_
        f_primer (str):
        r_primer (str): 
        blast_dir (Path): The BLAST bin path
        tmp_dir (Path): _description_
        threads (int): _description_
    """
    if not para_dict['tmp_dir'].exists():
        para_dict['tmp_dir'].mkdir()
    
    if not para_dict['out_path'].exists():
        para_dict['out_path'].mkdir()

    merged_seq_lst = merge_sanger_reads(para_dict['input_dir'], para_dict['extra_str'], para_dict['f_primer'], para_dict['r_primer'])
    query_seq_path = Path(para_dict['tmp_dir'] / 'all_sample_merged.fa')
    SeqIO.write(merged_seq_lst, query_seq_path,'fasta')
    logger.info('Merge Sanger sequencing reads. Done!')
    # Blast
    blast_cmd = NcbiblastnCommandline(
        cmd = para_dict['blast_dir'] / 'blastn',
        query = query_seq_path,
        db=para_dict['db_path'],
        evalue=1e-5,
        outfmt="6 qacc sacc qlen slen length qstart qend sstart send ssciname staxid pident evalue",
        max_hsps=1,
        max_target_seqs=20,
        num_threads=para_dict['threads'],
        out=para_dict['tmp_dir'] / 'blast.out'
    )
    blast_cmd()
    logger.info('BLAST against NCBI 16S rRNA gene database. Done!')
    # parse blast result
    tb_meta = pd.read_table(para_dict['meta_path'])
    tb_meta['queryid'] = tb_meta['row'] + tb_meta['column'].astype(str)
    tb_blast= pd.read_table(para_dict['tmp_dir'] / 'blast.out',
                  names=['queryid', 'subjectid', 'qlen', 'slen', 'alignlen',
                         'qstart', 'qend', 'sstart', 'send',
                         'species', 'taxid', 'identity', 'evalue'])
    tb_seq = pd.DataFrame([{'queryid': _seq.id, 'seq': _seq.seq} for _seq in merged_seq_lst])
    # merge
    tb_out = tb_seq.merge(tb_meta, how='left').merge(tb_blast, how='left')
    tb_out = tb_out[['row', 'column', 'plate', 'well', 'clone',
                     'subjectid', 'species', 'taxid', 'identity', 'alignlen', 'evalue',
                     'seq']]
    # Output
    # Output all result
    tb_out.to_csv(para_dict['out_path'] / 'all_hit.tsv', sep='\t', index=False)
    logger.info('Parse results. Done!')
    # Output the top hit for each sequence
    tb_dedup = tb_out.sort_values(by='evalue').drop_duplicates(subset=['row', 'column'])
    tb_dedup.to_csv(para_dict['out_path'] / 'top_tit.tsv',
                    sep='\t',
                    index=False)
    logger.info('Also output top hits. Done!')
    tb_dedup['seqid'] = tb_dedup['plate'].astype(str) + '_' + tb_dedup['well'].astype(str) + '_' + tb_dedup['clone'].astype(str)
    # construct the list of sequences
    seq_lst = []
    for _idx, _row in tb_dedup.iterrows():
        seq_lst.append(
            SeqRecord(id=_row['seqid'],
                    seq=Seq(_row['seq']),
                    description=''
                    )
        )
    keep_seq_lst = rm_dup_seqs(seq_lst, para_dict['tmp_dir'], para_dict['muscle_path'])
    keep_seq_id_lst = [_.id for _ in keep_seq_lst]
    logger.info('Output fasta sequence. Done!')
    # Also output non-redundant table
    tb_nonredundant = tb_dedup[tb_dedup['seqid'].isin(keep_seq_id_lst)]
    tb_nonredundant = tb_nonredundant[['row', 'column',	'plate', 'well', 'clone',
                                       'subjectid', 'species',	'taxid', 'identity',
                                       'alignlen', 'evalue', 'seqid','seq']]
    tb_nonredundant.to_csv(para_dict['out_path'] / 'Bacteria_ToIsolate.tsv',
                           sep='\t',
                           index=False)
    SeqIO.write(keep_seq_lst,
                para_dict['out_path'] / 'Bacteria_ToIsolate.fasta',
                'fasta')
    #shutil.rmtree(para_dict['tmp_dir'])

def parse_args():
    """Parse arguments
    """
    parser = argparse.ArgumentParser(description="Example with -1 and -2 options")
    parser.add_argument('-i', '--input', required=True, type=Path, dest='input_dir',
                        help='<file_path> Directory containing Sanger sequencing reads end with ".ab1".')
    parser.add_argument('-m', '--meta', required=True, type=Path, dest='meta_path',
                        help='<file_path> The meta info to match PCR plate and picked clones.')
    parser.add_argument('-d', '--db', required=True, type=Path, dest='db_path',
                        help='<file_path> The NCBI 16S rRNA gene database')
    parser.add_argument('-o', '--output', required=True, type=Path, dest='out_path',
                        help='<directory_path> The output directory path.')
    parser.add_argument('-f', '--fwr', type=str, dest='f_primer', default='27F',
                        help='<str> Forward primer name. Default: 27F')
    parser.add_argument('-r', '--rev', type=str, dest='r_primer', default='1492R',
                        help='<str> Reverse primer name. Default: 1492R')
    parser.add_argument('-e', '--extra', type=str, dest='extra_str', default='',
                        help='<str> Extra characters in file names that need to be removed')
    parser.add_argument('--blast', type=Path, dest='blast_dir', default=Path(''),
                        help='<file_path> The path to BLAST binary directory (if it is not in PATH)')
    parser.add_argument('--muscle', type=Path, dest='muscle_path', default=Path('muscle'),
                        help='<file_path> The path to MUSCLE5 (if it is not in PATH)')
    parser.add_argument('--threads', type=int, default=4,
                        help='<int> Threads for BLAST. Default: 4')
    parser.add_argument('--tmp', type=Path, dest='tmp_dir', default=Path('tmp'),
                        help='<int> Temperory directory path. Default: ./tmp')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    ARGS = parse_args()
    sanger_main(vars(ARGS))
