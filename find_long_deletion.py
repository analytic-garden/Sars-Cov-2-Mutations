#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
file_find_long_deleteion.py

@author: Bill Thompson
@license: GPL 3
@copyright: Sept. 8 2020
"""

import re
from consensus3 import ref_pos_to_alignment, read_alignment_file

def main():
    date = '2020_09_04'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    align_file =  base + 'sequences_no_dups_aln.fasta'
    deletion = '-' * 382
    ref_id = 'NC_045512.2'
    
    pos_map = ref_pos_to_alignment(align_file, ref_id)
    
    align = read_alignment_file(align_file)
    for id, record in align.items():
        for pos in [m.start() for m in re.finditer(deletion, str(record.seq))]:
        # pos = record.seq.find(deletion) 
            ref_pos = pos_map[pos]
            if ref_pos != -1 and ref_pos <= 29000:
                print(id, pos, pos_map[pos])

if __name__ == "__main__":
    main()

