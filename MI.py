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
    
    notes:
        This method is faster than using contingency_matrix from sklearn
    """
    nts = ['A', 'C', 'G', 'T']
    
    c1 = {'A': pseudo, 'C': pseudo, 'G': pseudo, 'T': pseudo}
    c2 = {'A': pseudo, 'C': pseudo, 'G': pseudo, 'T': pseudo}
    p = {'A': {'A': pseudo, 'C': pseudo, 'G': pseudo, 'T': pseudo},
          'C': {'A': pseudo, 'C': pseudo, 'G': pseudo, 'T': pseudo},
          'G': {'A': pseudo, 'C': pseudo, 'G': pseudo, 'T': pseudo},
          'T': {'A': pseudo, 'C': pseudo, 'G': pseudo, 'T': pseudo}}
    
    N = 0
    for i in range(len(col1)):
        if col1[i] in nts and col2[i] in nts:
            c1[col1[i]] += 1
            c2[col2[i]] += 1
            p[col1[i]][col2[i]] += 1
            N += 1
            
    mi = 0
    for nt1 in nts:
        count1 = c1[nt1]
        for nt2 in nts:
            count2 = c2[nt2]
            mi += (p[nt1][nt2]/N) * np.log2(N * p[nt1][nt2] / (count1 * count2))
            
    return mi
                 
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
    a dictionary
         'Position_1' - a list of genome positions
         'Position_2' - a list of genome positions
         'MI' - list of mutual information between Position_1 and Position_2
    
    requires:
    pseudo >= 0
    """
    var_cols = list(variant_cols.keys())
    
    mi_tab = {'Position_1': [], 'Position_2': [], 'MI': []}
    for i in range(len(var_cols)-1):
        for j in range(i+1, len(var_cols)):
            mi = MI(variant_cols[var_cols[i]][1],
                    variant_cols[var_cols[j]][1],
                    pseudo = pseudo)
            mi_tab['Position_1'].append(pos_map[var_cols[i]]+1)
            mi_tab['Position_2'].append(pos_map[var_cols[j]]+1)
            mi_tab['MI'].append(mi)
            
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
    mi_tab = MI_table(variant_cols, pos_map)
    df_mi = pd.DataFrame(mi_tab).sort_values(by = 'MI', ascending = False)
    df_mi.to_csv(mutual_info_csv, index=False)

if __name__ == "__main__":
    main()

