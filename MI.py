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
from consensus3 import MI_table, ref_pos_to_alignment, \
                        get_varying_columns, get_col_range

def main():
    date = '2020_10_22'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    align_file =  base + 'sequences_no_dups_aln.fasta'
    mutual_info_csv = base + 'MI_ncbi_no_dups.csv'

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

