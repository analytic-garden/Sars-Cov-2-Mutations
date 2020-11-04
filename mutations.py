#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mutations.py
    Find aligned columns with significant mutations 
    i.e. changes from reference seq.

@author: Bill Thompson
@license: GPL 3
@copyright: Oct. 19, 2020
"""

import pandas as pd
from Bio import AlignIO
from consensus3 import ref_pos_to_alignment, read_genbank_file, \
                        get_col_range, count_mutations
                        
def main():
    date = '2020_11_03'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    ## align_file =  base + 'sequences_no_dups_aln.fasta'
    align_file =  base + 'sequences_valid_dates_aln.fasta'
    genbank_file = base + 'ncbi_dataset/data/genomic.gbff'
    mutation_csv_base = base + 'sars_cov_2_ncbi_ncbi_mut_no_dups_'

    ref_id = 'NC_045512.2'

    # change these for different alignments and mutation level cutoffs
    min_col_quality = 0.90
    consensus_cutoff = 0.999
    mutation_csv_file = mutation_csv_base + str(consensus_cutoff * 100) + '.csv'

    # map alignment to the reference sequence
    pos_map = ref_pos_to_alignment(align_file, ref_id)

    # get the data
    align = AlignIO.read(align_file, 'fasta')
    genbank = read_genbank_file(genbank_file)
    
    start, end = get_col_range(align, min_col_quality)

    ref_seq = genbank[ref_id]  # the reference sequence record

    mutations = count_mutations(align, pos_map, ref_seq, 
                                start = start, end = end,
                                consensus_cutoff=consensus_cutoff)
    df2 = pd.DataFrame(mutations)
    df2 = df2.sort_values('consensus %')
    df2.to_csv(mutation_csv_file, index=False)

if __name__ == "__main__":
    main()

