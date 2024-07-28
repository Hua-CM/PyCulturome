# -*- coding: utf-8 -*-
# @Time    : 2024/7/14
# @Author  : Zhongyi Hua
# @File    : pick_wells.py
# @Usage   : pick wells based on asv counts
# @Note    : 
# @E-mail  : njbxhzy@hotmail.com
import argparse
from collections import defaultdict
from pathlib import Path

from ete3 import Tree, TreeStyle, NodeStyle, FaceContainer, TextFace, RectFace
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import seaborn as sns

# set plot parameters
mpl.rcParams['pdf.fonttype'] = 'TrueType'
mpl.rcParams['pdf.use14corefonts'] = True


def tidy_taxonomy(taxonomy_path: Path):
    """Tidy taxonomy annotation generated from SINTAX

    Args:
        taxonomy_path (Path): _description_

    Returns:
        _type_: _description_
    """
    tb_tax = pd.read_table(taxonomy_path, header=None)
    tb_tax_new = pd.DataFrame(tb_tax[0])
    tb_tax_new[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = tb_tax[3].str.split(',', expand=True).\
        fillna('Unclassified').\
            replace('[dpcofgs]:', '', regex=True)
    tb_tax_new.rename(columns={0: 'OTUID'}, inplace=True)
    return tb_tax_new


class SummaryLibrary:
    """_summary_
    """
    def __init__(self, df_count, control_pos, control_neg, threshold, taxonomy) -> None:
        self.df_count = df_count
        self.well_pos = control_pos
        self.well_neg = control_neg
        self.purity_thresh = threshold
        self.taxonomy = taxonomy
        # 
        self.df_positive = self.detect_pos_wells()
        self.df_purified = self.detect_purified_wells()
        self.df_asv_recommend = self.detect_recommend_wells()
        self.tree = self.construct_tree()

    def detect_pos_wells(self):
        """Detect wells with bacteria

        Args:
            df_count (pd.DataFrame): _description_
            control_neg (str): Negative control well id
            control_pos (str): Positive control well id
        """
        _df_count = self.df_count[['reads', 'well', 'plate']]
        _df_count = _df_count.groupby(['well', 'plate']).sum().reset_index()
        # If the number of reads is greater than 95% of the negative control wells, then the well is considered to contain bacteria
        control_reads_neg = _df_count.loc[_df_count['well'] == self.well_neg, 'reads'].to_list()

        threshold_reads_neg = np.quantile(control_reads_neg, 0.95)
        _df_count = _df_count.loc[_df_count['reads'] > threshold_reads_neg]
        _df_count = _df_count[_df_count['well'] != self.well_pos]
        return _df_count

    #Draw boxplot for positive and negative control
    def draw_control_boxplot(self):
        """Draw positive and negative control boxplot

        Args:
            df_data (pd.DataFrame): _description_
        """
        _df_count = self.df_count[['reads', 'well', 'plate']]
        _df_count = _df_count.groupby(['well', 'plate']).sum().reset_index()
        # If the number of reads is greater than 95% of the negative control wells,
        # then the well is considered to contain bacteria
        control_reads_pos = _df_count.loc[_df_count['well'] == self.well_pos, 'reads'].to_list()
        control_reads_neg = _df_count.loc[_df_count['well'] == self.well_neg, 'reads'].to_list()

        fig = plt.figure(figsize=(5,5), num=1)
        ax1 = fig.add_subplot()
        sns.boxplot([control_reads_pos, control_reads_neg],
                    ax=ax1)
        ax1.set(xlabel='Controls',
                ylabel='Count of reads',
                title='Positive and negative control',
                xticklabels=['Positive', 'Negative'])
        return fig

    def draw_reads_freq(self):
        """
        Returns:
            matplotlib.figure.Figure: The reads frequency plot
        """
        fig = plt.figure(figsize=(8, 6))
        ax = fig.subplots(1, 2)

        sns.histplot(self.df_positive['reads'],
                     bins=np.arange(min(self.df_positive['reads']),
                                    max(self.df_positive['reads']) + 500,
                                    500),
                ax=ax[0])
        ax[0].set(title='Frequency of reads per well',
                  xlabel='Count of reads',
                  ylabel='Count of wells')


        plate_data = self.df_positive.groupby('plate').sum('reads')['reads'].to_list()
        sns.histplot(
            plate_data,
            bins=np.arange(min(plate_data),
                           max(plate_data) + 10000,
                           500),
            ax=ax[1]
        )
        ax[1].set(title='Frequency of reads per plate',
                  xlabel='Count of reads',
                  ylabel='Count of plates')
        fig.set_facecolor('white')
        return fig


    def detect_purified_wells(self):
        """Detect purified wells for subsequent isolation

        Args:
            df_pos_wells (pd.DataFrame): _description_
            threshold (float): 0.95 means that a well is considered purified if 95% of 
            its reads belong to the same ASV.

        Returns:
            pd.DataFrame: A DataFrame with four columns: plate/well/OTUID/reads/purity/purified_well
        """
        self.df_positive['tmpwell'] = self.df_positive['plate'] + '_' + self.df_positive['well']
        self.df_count['tmpwell'] = self.df_count['plate'] + '_' + self.df_count['well']
        _df_count = self.df_count[self.df_count['tmpwell'].isin(self.df_positive['tmpwell'].to_list())]

        wells_dct = defaultdict(dict)
        for _idx, _row in _df_count.iterrows():
            wells_dct[f'{_row["plate"]}_{_row["well"]}'].update({_row['OTUID']: _row['reads']})
        # detect purified wells
        purified_well_lst = []
        for _well, _otu_dct in wells_dct.items():
            _max_key = max(_otu_dct, key=_otu_dct.get)
            _purity = (_otu_dct[_max_key] / sum(_otu_dct.values()))
            _puritified_well = bool(_purity > self.purity_thresh)
            purified_well_lst.append({'tmpwell': _well,
                                      'OTUID': _max_key,
                                      'purity': _purity,
                                      'purified_well': _puritified_well})
        df_purity = pd.DataFrame(purified_well_lst)
        
        df_purity = self.df_positive.merge(df_purity, how='left')
        df_purity.drop('tmpwell', axis=1, inplace=True)
        self.df_positive.drop('tmpwell', axis=1, inplace=True)
        self.df_count.drop('tmpwell', axis=1, inplace=True)
        return df_purity

    def draw_purity_freq(self):
        """Draw the frequency histogram of purified wells

        Returns:
            matplotlib.figure.Figure: the frequency histogram of purified wells.
        """
        def to_percent(_x, pos):
            return f'{_x * 100:.0f}%'

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot()
        sns.histplot(self.df_purified, x='purity', ax=ax)
        ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(to_percent))
        ax.set(title='Rarefication curve',
               xlabel='Purity(%)',
               ylabel='Well numbers')
        return fig

    def draw_rarefy_wells(self, n_bs=100):
        """Draw rarefication curve

        Args:
            n_bs (int, optional): _description_. Defaults to 100.

        Returns:
            _type_: _description_
        """
        _df_raw = self.df_purified[self.df_purified['purified_well']==True]
        
        _total_dct = {}
        for n_well in range(1, len(_df_raw), round(len(_df_raw) / 30)):
            n_well_lst = [
                len(set(_df_raw.sample(n_well)['OTUID'].to_list())) for i in range(n_bs)
            ]
            _total_dct[n_well] = n_well_lst
        _df_plot = pd.DataFrame(_total_dct)
        _df_plot = _df_plot.reset_index().melt(id_vars=['index'],
                                               var_name='numbers',
                                               value_name='ASVs')

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot()
        sns.boxplot(_df_plot, x='numbers', y='ASVs', ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, va='top')
        ax.set(title='Rarefication curve',
               xlabel='Well numbers',
               ylabel='ASV number')
        return fig


    def detect_recommend_wells(self):
        """Recommend wells for each ASV

        Returns:
            _type_: _description_
        """
        _df_raw = self.df_purified[self.df_purified['purified_well']==True]
        _df_res = _df_raw.groupby('OTUID').max(['purity', 'reads']).reset_index()
        return _df_res


    def construct_tree(self):
        """Construct 

        Returns:
            ete3.Tree: The taxonomic tree
        """
        tb_tree = self.taxonomy[self.taxonomy['OTUID'].isin(self.df_purified['OTUID'].to_list())]
        tax_tree = Tree()
        for _idx, _row in tb_tree.iterrows():
            _otu = _row[0]
            parent_node = tax_tree
            for _tax_name in _row[1:7]: # Only kingdom to genus
                _cur_node = tax_tree.search_nodes(name=_tax_name)
                if not _cur_node:
                    _cur_node = parent_node.add_child(name=_tax_name)
                else:
                    _cur_node = _cur_node[0]
                parent_node = _cur_node
            parent_node.add_child(name=_otu)
        return tax_tree


    def draw_tree(self, output_path: Path):
        # generate colors
        def generate_multiple_colors(num_colors):
            base_cmap = plt.get_cmap('tab20')  # Basic color palette
            colors = []
            for i in range(num_colors):
                colors.append(
                    mpl.colors.rgb2hex(
                        # Recycling colors and Set alpha to 0.3
                        base_cmap(i % base_cmap.N)[0:3] + (0.3,)
                        )
                    )
            return colors
        phylums = [_node.name for _node in self.tree.get_children()[0].get_children()]
        phylum_color_dct = dict(zip(phylums, generate_multiple_colors(len(phylums))))

        def set_node_style(node, color):
            style = NodeStyle()
            style["bgcolor"] = color
            style["size"] = 10
            style["vt_line_width"] = 2
            style["hz_line_width"] = 2
            node.set_style(style)

        def custom_layout(node):
            phylum_col = phylum_color_dct.get(node.name, '#FFFFFF')
            set_node_style(node, phylum_col)
            if node.is_leaf():
                setattr(node.faces, 'aligned', FaceContainer())
                face = TextFace(node.name, fsize=24)  # 这里 fsize 是字体大小，例如 14
                node.add_face(face, column=0, position='aligned')
        
        ts = TreeStyle()
        # draw legend
        t_legend = FaceContainer()
        for _phy, _col in phylum_color_dct.items():
            t_legend.add_face(RectFace(320, 160, '#FFFFFF', _col), column=0)
            t_legend.add_face(TextFace(_phy, fsize=60), column=1)
        ts.legend = t_legend
        ts.layout_fn = custom_layout
        ts.mode = 'c'
        ts.show_leaf_name = False
        self.tree.render(str(output_path), tree_style=ts)
        

def parse_args():
    """Parse CLI commands
    """
    parser = argparse.ArgumentParser(description="Assign plates and wells name to the sequence \
                                     id for subsequence analysis")
    parser.add_argument('-i', '--input', required=True, type=Path,
                        help='<file_path>  The asv count table generated by vsearch')
    parser.add_argument('-t', '--taxonomy', required=True, type=Path,
                        help='<file_path>  The asv taxonomy table generated by vsearch')
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help='<directory_path> The output directory path.')
    parser.add_argument('-p', '--positive', type=str, default='B12',
                        help='<str> Positive control well ID. Default: A12')
    parser.add_argument('-n', '--negative', type=str, default='A12',
                        help='<str> Negative control well ID. Default: B12')
    parser.add_argument('--threshold', type=float, default=0.95,
                        help='<float> Threshold for determining purified wells. Default: 0.95\
                        means that a well is considered purified if 95% of its reads belong to the same ASV.')
    args = parser.parse_args()
    return args


def summary_library(tb_abun_path: Path,
                    taxonomy_path: Path,
                    output_path: Path,
                    positive_well: str,
                    negative_well: str,
                    threshold: float
                    ):
    """Summary library features

    Args:
        tb_abun_path (Path): The usearch output file
        taxonomy_path (Path): The usearch sintax output file
        output_path (Path): The output directory
        positive_well (str):Positive control well ID
        negative_well (str): Negative control well ID
        threshold (float): The threshold for a 'purified' well
    """
    tb_count = pd.read_table(tb_abun_path)
    tb_count = tb_count.melt(id_vars=['#OTU ID'],
                             var_name='plate_well',
                             value_name='reads')
    tb_count[['well', 'plate']] = tb_count['plate_well'].str.split('_', expand=True)
    tb_count = tb_count[(tb_count['plate'] != 'Unclassified') & (tb_count['well'] != 'Unclassified')]
    tb_count.drop(columns=['plate_well'], inplace=True)
    tb_count.rename(columns={'#OTU ID': 'OTUID'}, inplace=True)
    tb_count[tb_count['reads'] != 0].to_csv(output_path / 'readscount_well.tsv', sep='\t', index=False)
    
    tb_tax = tidy_taxonomy(taxonomy_path)
    tb_tax.to_csv(output_path / 'taxonomy_8_asv.tsv', sep='\t', index=False)

    IsoLibrary = SummaryLibrary(tb_count, positive_well, negative_well, threshold, tb_tax)
    tb_wells_purified = IsoLibrary.df_purified.merge(IsoLibrary.taxonomy, how='left')

    # Table
    # # Purified wells
    tb_wells_purified.to_csv(output_path / 'purified_well.tsv', sep='\t', index=False)
    # # Recommend wells for each ASV
    IsoLibrary.df_asv_recommend.to_csv(output_path / 'recommend_well4asv.tsv', sep='\t', index=False)  

    # Plot
    fig1 = IsoLibrary.draw_control_boxplot()
    fig2 = IsoLibrary.draw_reads_freq()
    fig3 = IsoLibrary.draw_purity_freq()
    fig4 = IsoLibrary.draw_rarefy_wells()
    fig1.savefig(output_path / 'Positive_Negative_Control.pdf')
    fig2.savefig(output_path / 'Reads_frequency.pdf')
    fig3.savefig(output_path / 'Purity_frequency.pdf')
    fig4.savefig(output_path / 'Rarefication.pdf')
    # This function output to file directly because it relies on PyQt5
    IsoLibrary.draw_tree(output_path / 'Tree.pdf')

def main():
    """Main interface
    """
    args = parse_args()
    if not args.output.exists():
        args.output.mkdir()
    summary_library(args.input,
                    args.taxonomy,
                    args.output,
                    args.positive,
                    args.negative,
                    args.threshold)

if __name__ == '__main__':
    main()
