#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
variation.py
    find locataions with significant variation

@author: Bill Thompson
@license: GPL 3
@copyright: Oct. 19, 2020
"""

from Bio import AlignIO
from consensus3 import read_genbank_file, init_dataframe, get_varying_columns, \
                        get_col_range, count_variation_from_ref, \
                        ref_pos_to_alignment
                        
def main():
    date = '2020_11_11'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    ## align_file =  base + 'sequences_no_dups_aln.fasta'
    align_file =  base + 'sequences_valid_dates_aln.fasta'
    csv_file = base + 'sars_cov_2_variation_ncbi_valid_dates.csv'
    genbank_file = base + 'ncbi_dataset/data/genomic.gbff'

    ref_id = 'NC_045512.2'
    min_col_quality = 0.90
    
    # get the data
    align = AlignIO.read(align_file, 'fasta')
    genbank = read_genbank_file(genbank_file)

    # change these for different alignments and mutation level cutoffs
    
    # map alignment to the reference sequence
    pos_map = ref_pos_to_alignment(align_file, ref_id)
    
    start, end = get_col_range(align, min_col_quality)

    df = init_dataframe(genbank, align)
        
    # get varying columns
    variant_cols = get_varying_columns(align, 
                                       start = start, end = end)

    diffs = count_variation_from_ref(align, align_file, ref_id, start, end)
    df['Diffs from Ref'] = diffs
    
    # wite the data
    for col in variant_cols:
        df['Pos ' + str(pos_map[col]+1)] = variant_cols[col][1]

    df = df.sort_values('Collection Date')
    df.to_csv(csv_file, index=False)

if __name__ == "__main__":
    main()

