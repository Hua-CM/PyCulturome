# -*- coding: utf-8 -*-
# @File    :   validate_asv.py
# @Time    :   2024/09/18 21:23:27
# @Author  :   Zhongyi Hua 
# @Usage   :   Validate whether the isolated strain is identical to the sequenced ASV
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com
import logging
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

from .utils import pairwise_align, runCommand

logger = logging.getLogger(__name__)
pd.options.display.float_format = '{:,.0f}'.format


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
        sanger_path (Path): The summarized Sanger sequencing result path. MUST WITH *seq* column
        ngs_path (Path): The NGS analysis result: purified_well.tsv.
        asv_path (Path): The *.rep.fa.
        out_path (Path): The result output path.
        blast_dir(Path): BLAST software directory
        tmp_dir  (Path): temporary directory
        threads  (int) : Number of threads
    """
    tb_sanger = pd.read_table(para_dict['sanger_path'])
    tb_ngs = pd.read_table(para_dict['ngs_path'])
    asv_seq_dct = SeqIO.to_dict(SeqIO.parse(para_dict['asv_path'], 'fasta'))
    well_dct = {}
    for _idx, _row in tb_ngs.iterrows():
        well_dct.setdefault(_row['plate'] + _row['well'], _row['OTUID'])
    check_res_lst = []
    inconsistent_seq_lst = []
    for _idx, _row in tb_sanger.iterrows():
        raw_well = _row['plate'] + _row['well']
        try:
            seq_sanger = _row['seq']
            seq_ngs = asv_seq_dct.get(well_dct.get(raw_well))
            alignment = pairwise_align(seq_sanger, seq_ngs)
            consistent = check_identical(alignment[0], alignment[1])
            check_res_lst.append({'plate': _row['plate'],
                                  'well' : _row['well'],
                                  'clone': _row['clone'],
                                  'consistency': consistent})
            if consistent == 0:
                inconsistent_seq_lst.append(
                    SeqRecord(
                        id = '_'.join([_row['plate'], _row['well'], str(_row['clone'])]),
                        seq=Seq(seq_sanger),
                        description= '',
                        name=''
                    )
                )
        except:
            logger.info(f'{raw_well} not in NGS result, please check')
    tb_validate = pd.DataFrame(check_res_lst)
    tb_sanger = tb_sanger.merge(tb_ngs[['plate', 'well', 'OTUID']], how='left').merge(tb_validate)

    if inconsistent_seq_lst:
    #  check inconsistent squences 
    # #Extract inconsistent sequences
        if not para_dict['tmp_dir'].exists():
            para_dict['tmp_dir'].mkdir()
        run_makeblastdb(para_dict['asv_path'],
                        para_dict['tmp_dir'] / 'asv_db',
                        str(para_dict['blast_dir'] / 'makeblastdb'))
        query_seq_path = para_dict['tmp_dir'] / 'query.fasta'
        SeqIO.write(inconsistent_seq_lst, query_seq_path, 'fasta')
        blastn_para_dct = {
            'bin': str(para_dict['blast_dir'] / 'blastn'),
            'evalue': 1e-5,
            'outfmt': '6 qacc sacc slen length pident evalue',
            'max_hsps': 1,
            'max_target_seqs': 1,
            'num_threads': str(para_dict['threads']),
        }
        run_blastn(query_seq_path,
                para_dict['tmp_dir'] / 'asv_db',
                para_dict['tmp_dir'] / 'asv_res.tsv',
                blastn_para_dct)
        # parse Sequence query against ASV result
        tb_res_asv= pd.read_table(para_dict['tmp_dir'] / 'asv_res.tsv',
                                names=['queryid', 'subjectid', 'slen', 'alignlen',
                                        'identity', 'evalue'])
        inconsistent_res_lst = []
        use_ncbi_lst = []
        for _idx, _row in tb_res_asv.iterrows():
            if (_row['slen'] == _row['alignlen']) and (_row['identity'] == 100):
                inconsistent_res_lst.append({'seqid': _row['queryid'],
                                            'TrueOTUID': _row['subjectid']})
            else:
                use_ncbi_lst.append(_row['queryid'])
        # Use sanger NCBI result directly.
        tb_inconsistent = pd.DataFrame(inconsistent_res_lst)
        tb_sanger = tb_sanger.merge(tb_inconsistent, how='left')
        if use_ncbi_lst:
            tb_sanger.loc[tb_sanger['seqid'].isin(use_ncbi_lst), 'TrueOTUID'] = tb_sanger.loc[tb_sanger['seqid'].isin(use_ncbi_lst), 'subjectid']
        tb_sanger.loc[tb_sanger['TrueOTUID'].isna(), 'TrueOTUID'] = tb_sanger.loc[tb_sanger['TrueOTUID'].isna(), 'OTUID']
    else:
        tb_sanger['TrueOTUID'] = tb_sanger['OTUID']

    # move seq column to the last
    columns = list(tb_sanger.columns)
    columns.append(columns.pop(columns.index('seq')))
    tb_sanger = tb_sanger[columns]
    tb_sanger.to_csv(para_dict['out_path'], sep='\t', index=False)
