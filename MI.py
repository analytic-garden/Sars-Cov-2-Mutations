#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MI.py
    Calculate mutual information among aligned columns

@author: Bill Thompson
@license: GPL 3
@copyright: Oct. 19, 2020
"""

from Bio import AlignIO
import pandas as pd
from sklearn import metrics
import numpy as np
from consensus3 import ref_pos_to_alignment, \
                        get_varying_columns, get_col_range


def MI(col1, col2, pseudo = 0.05):
    """
    MI - calculate the mutual information between two columns of a 
         multisequence alignment
    
    arguments:
    col1, col2 - two columns from the alignment
    pseudo - pseudocounts
    
    returns:
    the mutual information in bits
    
    requires:
    len(col1) == len(col2)
    pseudo >= 0
    np.log2
    """
    return metrics.mutual_info_score(col1, col2)/np.log(2)

def MI_table(variant_cols, pos_map, pseudo = 0.05):
    """
    MI_table - generate a table of mutual information values from all paires of
               columns in a  dictionary of aligned nucleotides
    
    arguments:
    variant_cols - a dictionary of columns. 
                   key = aligned column number, value = column data
    pos_map - a dictionary, key = aligned column position, value = reference coulum  position
    pseudo - pseudocounts
    
    returns:
    a dictionary -  sorted by value
         key = string "col1_number, col2_number", 
         value = mutual information between columns
    
    requires:
    pseudo >= 0
    """
    var_cols = list(variant_cols.keys())
    mi_tab = dict()
    for i in range(len(var_cols)-1):
        for j in range(i+1, len(var_cols)):
            k = str(pos_map[var_cols[i]]+1) + ', ' + str(pos_map[var_cols[j]]+1)
            mi_tab[k] = MI(variant_cols[var_cols[i]][1],
                           variant_cols[var_cols[j]][1],
                           pseudo = pseudo)
    mi_tab = {k: v for k, v in sorted(mi_tab.items(), key=lambda item: item[1], reverse=True)}

    return mi_tab

def main():
    date = '2020_11_11'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    ## align_file =  base + 'sequences_no_dups_aln.fasta'
    align_file =  base + 'sequences_valid_dates_aln.fasta'
    mutual_info_csv = base + 'MI_ncbi_valid_dates.csv'

    min_col_quality = 0.90
    consensus_cutoff = 0.98
    ref_id = 'NC_045512.2'

    # change these for different alignments and mutation level cutoffs
    mi_cutoff = 0.5
    
    # get the data
    align = AlignIO.read(align_file, 'fasta')
 
    # map alignment to the reference sequence
    pos_map = ref_pos_to_alignment(align_file, ref_id)
    
    start, end = get_col_range(align, min_col_quality)

   # get varying columns
    variant_cols = get_varying_columns(align, 
                                       consensus_cutoff = consensus_cutoff,
                                       start = start, end = end)

    # mutual information
    print('Mutual Infomation')
    mi_tab = MI_table(variant_cols, pos_map)
    mi_dict = {'Position_1': [], 'Position_2': [], 'MI': []}
    for k,mi in mi_tab.items():
        p1, p2 = k.split(',')
        mi_dict['Position_1'].append(int(p1))
        mi_dict['Position_2'].append(int(p2))
        mi_dict['MI'].append(mi)
        if mi >= mi_cutoff:
            print('MI(', k,  ') =', mi)
    df_mi = pd.DataFrame(mi_dict)
    df_mi.to_csv(mutual_info_csv, index=False)


if __name__ == "__main__":
    main()

