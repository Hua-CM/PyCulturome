# -*- coding: utf-8 -*-
# @File    :   GUI.py
# @Time    :   2024/09/02 22:34:21
# @Author  :   Zhongyi Hua
# @Usage   :   Make a GUI for Yuxuan Zhang
# @Note    :   
# @E-mail  :   njbxhzy@hotmail.com
import logging
import sys
from pathlib import Path
import traceback

from Bio.Align import fasta # For pyinstaller wrapper

import PySimpleGUI as sg
import tkinter as tk

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from PyCulturome.downstream.sanger import sanger_main
from PyCulturome.downstream.update_database import update_db_main
from PyCulturome.downstream.rename import rename_main
from PyCulturome.downstream.validate_asv import validate_asv
from PyCulturome.downstream.check_sanger import check_sanger_main


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = Path(sys._MEIPASS)
    except Exception:
        base_path = Path(".").resolve()

    return base_path / relative_path

image_path = resource_path("icon10.png")
ico_path = resource_path("meaw.ico")


def rename_seq_tab():
    """The tab for renaming sequences for long-term preservation
    """
    layout = [[sg.Frame('Data setting', [
            [sg.T('The input fasta path:',  size=25), sg.I('a fasta file', key='-RENAME INPUT-', size=35), sg.FileBrowse(target='-RENAME INPUT-')],
            [sg.T('The output fasta path:',  size=25), sg.I('a fasta file', key='-RENAME OUTPUT-', size=35), sg.FileBrowse(target='-RENAME OUTPUT-')],
            [sg.T('The old and new sequence id:', size=25), sg.I('Must has "oldid" and "newid" columns', key='-RENAME META-', size=35), sg.FileBrowse(target='-RENAME META-')],
             ])]]
    return sg.Tab('Rename', layout=layout, expand_x = True)

def update_db_tab():
    """The tab for updating existing bacteria database
    """
    row1 =sg.Frame('Data setting', [
            [sg.T('The existing database path:',  size=25), sg.I('a fasta file', key='-UPDATE DB-', size=35), sg.FileBrowse(target='-UPDATE DB-')],
            [sg.T('The sequences to add:',  size=25), sg.I('a fasta file', key='-UPDATE INPUT-', size=35), sg.FileBrowse(target='-UPDATE INPUT-')],
            [sg.T('The meta for sequences to add:', size=25), sg.I('Must has a column named sequence id', key='-UPDATE META-', size=35), sg.FileBrowse(target='-UPDATE META-')],
            [sg.T('Final sequences to be added:', size=25), sg.I(key='-UPDATE OUTPUT-', size=35), sg.FolderBrowse(target='-UPDATE OUTPUT-')]
             ])
    row2 = sg.Frame('Global setting', [
            [sg.T('Temporary directory path (Optional):', size=30, pad=(0,0)),
             sg.I('tmp', key='-UPDATE TMP_DIR-', size=30),
             sg.FolderBrowse(target='-UPDATE TMP_DIR-')],
            [sg.T('MUSCLE binary path (Optional):', size=30, pad=(0,0)), sg.I(key='-UPDATE MUSCLE_BIN-', size=30), sg.FolderBrowse(target='-UPDATE MUSCLE_BIN-')],
        ])
    main_col = [
        [row1],
        [row2]
    ]
    return sg.Tab('Update', layout=main_col, expand_x = True)


def sanger_tab():
    """The tab for analyzing Sanger sequencing result
    """
    row1 =sg.Frame('Data setting', [
            [sg.T('Sanger reads directory:',  size=20), sg.I(key='-SANGER INPUT-', size=40), sg.FolderBrowse(target='-SANGER INPUT-')],
            [sg.CB('Merge Sanger reads bebased on sample', key='-SANGER MERGE-', default=True, size=40)],
            [sg.T('Forward primer name:', size=15), sg.I('27F', key='-PRIMER FORWARD-', size=10),
             sg.T('Reverse primer name:', size=15), sg.I('1492R', key='-PRIMER REVERSE-', size=10)],
            [sg.T('Extra characters in file names that need to be removed',  size=40), sg.I(key='-SANGER EXTRA-', size=25)],
            [sg.T('Meta table path:', size=20), sg.I(key='-SANGER META-', size=40), sg.FileBrowse(target='-SANGER META-')],
            [sg.CB('Use "sample" as id', key='-SANGER ID1-', enable_events=True, default=False, size=30),
             sg.CB('Use "row" and "column" as id', key='-SANGER ID2-', enable_events=True, default=True, size=30)],
            [sg.T('Output directory:', size=20), sg.I(key='-SANGER OUTPUT-', size=40), sg.FolderBrowse(target='-SANGER OUTPUT-')],
            [sg.T('NCBI 16S Database path:', size=20), sg.I(key='-SANGER DATABASE-', size=40), sg.FileBrowse(target='-SANGER DATABASE-')]
             ])
    row2 = sg.Frame('Global setting', [
            [sg.T('Temporary directory path (Optional):', size=30, pad=(0,0)),
             sg.I('tmp', key='-SANGER TMP_DIR-', size=30),
             sg.FolderBrowse(target='-SANGER TMP_DIR-')],
            [sg.T('BLAST directory (Optional):', size=30, pad=(0,0)), sg.I(key='-SANGER BLAST_DIR-', size=30), sg.FolderBrowse(target='-SANGER BLAST_DIR-')],
            [sg.T('MUSCLE binary path (Optional):', size=30, pad=(0,0)), sg.I(key='-SANGER MUSCLE_BIN-', size=30), sg.FileBrowse(target='-SANGER MUSCLE_BIN-')],
            [sg.T('Threads (Optional):', size=30, pad=(0,0)), sg.I('4', key='-SANGER THREADS-')]
        ])
    main_col = [
        [row1],
        [row2]
    ]
    return sg.Tab('Sanger', layout=main_col, expand_x = True)


def validate_tab():
    """The tab for updating existing bacteria database
    """
    row1 =sg.Frame('Data setting', [
            [sg.T('The Sanger sequencing result path:',  size=25), sg.I('The summarised sanger table', key='-VALI SANGER-', size=35), sg.FileBrowse(target='-VALI SANGER-')],
            [sg.T('The NGS plate-well-ASV info:',  size=25), sg.I('a TSV file with plate/well/OTUID', key='-VALI NGS-', size=35), sg.FileBrowse(target='-VALI NGS-')],
            [sg.T('The NGS sequence path:', size=25), sg.I('The ASV sequence fasta path', key='-VALI SEQ-', size=35), sg.FileBrowse(target='-VALI SEQ-')],
            [sg.T('Output table path:', size=25), sg.I(key='-VALI OUTPUT-', size=35), sg.FolderBrowse(target='-VALI OUTPUT-')]
             ])
    row2 = sg.Frame('Global setting', [
            [sg.T('Temporary directory path (Optional):', size=30, pad=(0,0)),
             sg.I('tmp', key='-VALI TMP_DIR-', size=30),
             sg.FolderBrowse(target='-VALI TMP_DIR-')],
            [sg.T('BLAST directory (Optional):', size=30, pad=(0,0)), sg.I(key='-VALI BLAST_DIR-', size=30), sg.FolderBrowse(target='-VALI BLAST_DIR-')],
            [sg.T('Threads (Optional):', size=30, pad=(0,0)), sg.I('4', key='-VALI THREADS-')]
        ])
    main_col = [
        [row1],
        [row2]
    ]
    return sg.Tab('Validate', layout=main_col, expand_x = True)


def compare_tab():
    row1 =sg.Frame('Data setting', [
            [sg.T('The old squence table path:',  size=25), sg.I('With seqid and seq column', key='-COMP OLD-', size=35), sg.FileBrowse(target='-COMP OLD-')],
            [sg.CB('Old sequence use reverse complement sequence', key='-COMP OLDRP-', default=False, size=40)],
            [sg.T('The new squence table path:',  size=25), sg.I('With seqid and seq column', key='-COMP NEW-', size=35), sg.FileBrowse(target='-COMP NEW-')],
            [sg.CB('New sequence use reverse complement sequence', key='-COMP NEWRP-', default=False, size=40)],
            [sg.T('Output directory path:', size=25), sg.I(key='-COMP OUTPUT-', size=35), sg.FolderBrowse(target='-COMP OUTPUT-')]
             ])
    main_col = [
        [row1]
    ]
    return sg.Tab('Compare', layout=main_col, expand_x = True)

def make_window():
    """
    Make the main window
    """
    sg.theme('SystemDefaultForReal')
    tab1 = sanger_tab()
    tab2 = validate_tab()
    tab3 = update_db_tab()
    tab4 = rename_seq_tab()
    tab5 = compare_tab()

    layout = [
        [sg.TabGroup([[tab1, tab2, tab3, tab4, tab5]], key='-TASK-')],
        [sg.ML('',
               size=(60,8),
               k='-OUT-',
               write_only=True,
               expand_x=True,
               reroute_stdout=True, #True
               reroute_stderr=True, #True
               echo_stdout_stderr=True, #True
               reroute_cprint=True,
               auto_refresh=True,
               autoscroll=True,
               rstrip=False)],
        [sg.Button('Run', key='run'), sg.Exit()]
        ]
    window = sg.Window('PyCulturome-Sanger Sequencing', layout, default_element_size=(30,1), icon=ico_path)
    return window


def main():
    """
    The main interface for event loop
    """
    window = make_window()
    
    logger = logging.getLogger('PyCulturome')
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    while True:
        try:
            event, values = window.read()
            if event in (sg.WINDOW_CLOSED, 'Exit'):
                break
            if event.endswith('MULTIIN-'):
                # for multiline box input. use "-* MULTIIN-" for invisible input box.
                multi_key = ' '.join(event.split(' ')[:-1]) + '-'
                window[multi_key].Update('\n'.join(values[event].split(';')) + '\n', append=True)
            if event in ('-SANGER ID1-', '-SANGER ID2-'):
                for key in ('-SANGER ID1-', '-SANGER ID2-'):
                    if key != event:
                        window[key].update(False)
            if event == 'run':
                window['-OUT-'].update('')
                if values['-TASK-'] == 'Sanger':
                    # PySimpleGUI need convert variable type manually
                    para_dict = {
                        'is_merge': values['-SANGER MERGE-'],
                        'input_dir': Path(values['-SANGER INPUT-']),
                        'meta_path': Path(values['-SANGER META-']),
                        'extra_str': str(values['-SANGER EXTRA-']),
                        'f_primer': str(values['-PRIMER FORWARD-']),
                        'r_primer': str(values['-PRIMER REVERSE-']),
                        'db_path': Path(values['-SANGER DATABASE-']),
                        'out_path': Path(values['-SANGER OUTPUT-']),
                        'tmp_dir': Path(values['-SANGER TMP_DIR-']),
                        'blast_dir': Path(values['-SANGER BLAST_DIR-']),
                        'muscle_path': Path(values['-SANGER MUSCLE_BIN-']), #
                        'threads': int(values['-SANGER THREADS-'])
                    }
                    if values['-SANGER ID1-']: # use sample as ID
                        para_dict['keep_raw'] = True
                    else: # use row and column as ID
                        para_dict['keep_raw'] = False
                    sanger_main(para_dict)
                if values['-TASK-'] == 'Update':
                    para_dict = {
                        'in_path': Path(values['-UPDATE INPUT-']),
                        'out_path': Path(values['-UPDATE OUTPUT-']),
                        'db_path': Path(values['-UPDATE DB-']),
                        'tmp_dir': Path(values['-UPDATE TMP_DIR-']),
                        'bin_path': Path(values['-UPDATE MUSCLE_BIN-']),
                        'meta_path': Path(values['-UPDATE META-'])
                    }
                    update_db_main(para_dict)
                if values['-TASK-'] == 'Rename':
                    para_dict = {
                        'in_path': Path(values['-RENAME INPUT-']),
                        'out_path': Path(values['-RENAME OUTPUT-']),
                        'meta_path': Path(values['-RENAME META-'])
                    }
                    rename_main(para_dict)
                if values['-TASK-'] == 'Validate':
                    para_dict = {
                        'sanger_path': Path(values['-VALI SANGER-']),
                        'ngs_path': Path(values['-VALI NGS-']),
                        'asv_path': Path(values['-VALI SEQ-']),
                        'out_path': Path(values['-VALI OUTPUT-']),
                        'blast_dir': Path(values['-VALI BLAST_DIR-']),
                        'tmp_dir': Path(values['-VALI TMP_DIR-']),
                        'threads': Path(values['-VALI THREADS-'])
                    }
                    validate_asv(para_dict)
                if values['-TASK-'] == 'Compare':
                    para_dict = {
                        'old_path': Path(values['-COMP OLD-']),
                        'new_path': Path(values['-COMP NEW-']),
                        'old_rp': values['-COMP OLDRP-'],
                        'new_rp': values['-COMP NEWRP-'],
                        'out_path': Path(values['-COMP OUTPUT-'])
                    }
                    check_sanger_main(para_dict)
                img = tk.PhotoImage(file=image_path)
                print('Successful. Meow~')
                window['-OUT-'].widget.image_create(tk.INSERT, image=img)
        except Exception as e:
            exc_type, exc_value, exc_tb = sys.exc_info()
            # 格式化异常信息
            formatted_exception = ''.join(traceback.format_exception(exc_type, exc_value, exc_tb))
            print(formatted_exception)
            continue
    window.close()

if __name__ == '__main__':
    main()
