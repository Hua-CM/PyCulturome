# -*- coding: utf-8 -*-
# @File    :   downstream_pipe.py
# @Time    :   2024/09/10 21:02:58
# @Author  :   Zhongyi Hua 
# @Usage   :   
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com
import argparse
import logging
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from PyCulturome.downstream.sanger import sanger_main
from PyCulturome.downstream.update_database import update_db_main
from PyCulturome.downstream.rename import rename_main
from PyCulturome.downstream.validate_asv import validate_asv
from PyCulturome.downstream.check_sanger import check_sanger_main


class CustomFormatter(argparse.HelpFormatter):
    """Reduce the redundant metavar
    """
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string


def parse_args():
    """Parse arguments
    """
    parser = argparse.ArgumentParser(
        prog='PyCulturome',
        description='This script was for High-throughout sequencing ',
        formatter_class=CustomFormatter)
    sub_parser = parser.add_subparsers(title='Available', dest='command')
    sub_parser.required = True

    sanger_parser = sub_parser.add_parser(
        'sanger', help='summary all Sanger sequencing results in a directory')
    sanger_parser.add_argument('-i', '--input', required=True, type=Path, dest='input_dir',
                        help='<file_path> Directory containing Sanger sequencing reads end with ".ab1".')
    sanger_parser.add_argument('-m', '--meta', required=True, type=Path, dest='meta_path',
                        help='<file_path> The meta info to match PCR plate and picked clones.')
    sanger_parser.add_argument('-d', '--db', required=True, type=Path, dest='db_path',
                        help='<file_path> The NCBI 16S rRNA gene database')
    sanger_parser.add_argument('-o', '--output', required=True, type=Path, dest='out_path',
                        help='<directory_path> The output directory path.')
    sanger_parser.add_argument('--no_merge', type=bool, dest='is_merge', action='store_false', default=True,
                        help='Not merge Sanger reads based on sample name')
    sanger_parser.add_argument('-f', '--fwr', type=str, dest='f_primer', default='27F',
                        help='<str> Forward primer name. Default: 27F')
    sanger_parser.add_argument('-r', '--rev', type=str, dest='r_primer', default='1492R',
                        help='<str> Reverse primer name. Default: 1492R')
    sanger_parser.add_argument('-e', '--extra', type=str, dest='extra_str', default='',
                        help='<str> Extra characters in file names that need to be removed')
    sanger_parser.add_argument('--blast', type=Path, dest='blast_dir', default=Path(''),
                        help='<file_path> The path to BLAST binary directory (if it is not in PATH)')
    sanger_parser.add_argument('--muscle', type=Path, dest='muscle_path', default=Path('muscle'),
                        help='<file_path> The path to MUSCLE5 (if it is not in PATH)')
    sanger_parser.add_argument('--threads', type=int, default=4,
                        help='<int> Threads for BLAST. Default: 4')
    sanger_parser.add_argument('--tmp', type=Path, dest='tmp_dir', default=Path('tmp'),
                        help='<int> Temperory directory path. Default: ./tmp')
    
    validate_parser = sub_parser.add_parser(
        'validate', help='Confirm whether the sequence of the isolated strain is consistent with the results of NGS.')
    validate_parser.add_argument('-i', '--sanger', required=True, type=Path, dest='sanger_path',
                        help='<file_path> The Sanger sequencing result path')
    validate_parser.add_argument('-n', '--ngs', required=True, type=Path, dest='ngs_path',
                        help='<file_path> The NGS plate-well-ASV info')
    validate_parser.add_argument('-s', '--seq', required=True, type=Path, dest='asv_path',
                        help='<file_path> The NGS sequence fasta path')
    validate_parser.add_argument('-o', '--out', required=True, type=Path, dest='out_path',
                        help='<file_path> The output table path')
    validate_parser.add_argument('--threads', type=int, default=4,
                        help='<int> Threads for BLAST. Default: 4')
    validate_parser.add_argument('--tmp', type=Path, dest='tmp_dir', default=Path('tmp'),
                        help='<int> Temperory directory path. Default: ./tmp')
    validate_parser.add_argument('--blast', type=Path, dest='blast_dir', default=Path(''),
                        help='<file_path> The path to BLAST binary directory (if it is not in PATH)')
    
    update_parser = sub_parser.add_parser(
        'update', help='Update the exsiting bacteria culture collection')
    update_parser.add_argument('-d', '--database', required=True, type=Path, dest='db_path',
                               help='<file_path> The existing dulture collection path. A fasta file.')
    update_parser.add_argument('-i', '--input', required=True, type=Path, dest='in_path',
                            help='<file_path> The sequences to add. A fasta file')
    update_parser.add_argument('-o', '--output', required=True, type=Path, dest='out_path',
                            help='<file_path> Final Non-redundant sequences to add. A fasta file')
    update_parser.add_argument('--meta', type=Path, dest='bin_path',
                        help='<file_path> meta information for sequences to add (If any). MUST WITH A "seqid" column.')
    update_parser.add_argument('-t', '--tmp', type=Path, dest='tmp_dir', default=Path('tmp'),
                            help='<file_path> The sequences to add. A fasta file')
    update_parser.add_argument('--muscle', type=Path, dest='bin_path', default=Path('muscle'),
                        help='<file_path> The path to MUSCLE5 (if it is not in PATH)')


    rename_parser = sub_parser.add_parser(
        'rename', help='Rename the sequence ID to the long-term preservation ID.')
    rename_parser.add_argument('-i', '--input', required=True, type=Path, dest='in_path',
                        help='<file_path> The raw sequence fasta path.')
    rename_parser.add_argument('-m', '--meta', required=True, type=Path, dest='meta_path',
                        help='<file_path> A TSV file contanins "oldid" and "newid".')
    rename_parser.add_argument('-o', '--output', required=True, type=Path, dest='out_path',
                        help='<file_path>  The converted sequence fasta path.')    
    
    compare_parser = sub_parser.add_parser(
       'compare', help='Compare the results of the two Sanger sequencing experiments.' 
    )
    compare_parser.add_argument('--old', required=True, type=Path, dest='old_path',
                        help='<file_path>  The old Sanger sequencing result table path. MUST WITH "seqid" and "seq" column')
    compare_parser.add_argument('--new', required=True, type=Path, dest='new_path',
                        help='<file_path> The new Sanger sequencing result table path. MUST WITH "seqid" and "seq" column')
    compare_parser.add_argument('--r1', type=bool, dest='old_rp', action='store_true', default=False,
                        help='Old sequence use reverse complement sequence')
    compare_parser.add_argument('--r2',  type=bool, dest='new_rp', action='store_true', default=False,
                        help='New sequence use reverse complement sequence')   
    compare_parser.add_argument('-o', '--output', required=True, type=Path, dest='out_path',
                        help='<file_path>  The output directory path.')  

    args = parser.parse_args()
    return args


def main(para_dct):
    logger = logging.getLogger('PyCulturome')
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    func_dct = {'sanger': sanger_main,
                'validate': validate_asv,
                'update': update_db_main,
                'rename': rename_main}
    func = func_dct.get(para_dct['command'])
    func(para_dct)


if __name__ == '__main__':
    ARGS = parse_args()
    main(vars(ARGS))
