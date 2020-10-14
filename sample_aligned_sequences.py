#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sample_aligned_sequences.py 
@description: sample a series if sequences from a FASTA alignment

@author: Bill Thompson
@license: GPL 3
@copyright: August 17, 2020
"""

import random
from consensus3 import read_alignment_file
from Bio import Align
from Bio import AlignIO

def main():
    date = '2020_10_07'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    in_file = base + 'sequences_filtered_aln2_red.fasta'
    out_file = base + 'trees/sequences_filtered_aln2_samp1.fasta'
    
    num_samples = 1000
    ref_id = 'NC_045512.2'
    
    aln = read_alignment_file(in_file)
    seq_ids = random.sample(list(aln.keys()), k=num_samples)
    if not ref_id in seq_ids:
        seq_ids[0] = ref_id
        
    alignment = Align.MultipleSeqAlignment([])
    for id in seq_ids:
       alignment.append(aln[id]) 
                
    AlignIO.write(alignment, open(out_file, 'w'), 'fasta')

if __name__ == "__main__":
    main()




