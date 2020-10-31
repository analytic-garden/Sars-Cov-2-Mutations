#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_alignment_cols.py
    create a new alignment from aligned column regiom

@author: Bill Thompson
@license: GPL 3
@copyright: Oct. 31, 2020
"""

from Bio import AlignIO

def main():
    date = '2020_10_27'
    base = '/mnt/g/Covid-19/' + date + '/' 

    align_in_file = base + 'sequences_valid_dates_aln.fasta'
    align_out_file = base + 'sequences_cols_aln.fasta'
    
    align_in = AlignIO.read(align_in_file, 'fasta')
    AlignIO.write(align_in[:, 14881:14884], 
                  align_out_file, 
                  'fasta')

if __name__ == "__main__":
    main()
