#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
trim_alignment_py - remove low quality regions from either end of alignment

@author: Bill Thompson
@license: GPL 3
@copyright: Aug. 14, 2020
"""

from Bio import AlignIO
from consensus3 import get_col_range

def main():
    date = '2020_10_07'
    base = '/mnt/g/Covid-19/' + date + '/' 
     
    align_file = base + 'sequences_filtered_aln.fasta'
    out_file = base + 'sequences_filtered_aln_red.fasta'
    
    align = AlignIO.read(open(align_file), 'fasta')
    start, end = get_col_range(align, min_quality = 0.90)
    
    align2 = align[:, start:(end+1)]
    AlignIO.write(align2, open(out_file, 'w'), 'fasta')

if __name__ == "__main__":
    main()

