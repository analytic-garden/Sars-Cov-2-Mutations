"""
filter_seqs.py - remove duplicate sequences and 
                  add date and location to FASTA header
author: Bill Thompson
license: GPL 3
copyright: Apr. 30, 2020
"""

import os
from Bio import SeqIO
import pandas as pd

"""
read_fasta_file - read a collection sequences in FASTA format

arguments:
filename - the name of the FASTA file

returns:
a dictionary. The keys are sequence IDs. The values are BioPython SeqRecords.
"""
def read_fasta_file(filename):
    seqs = SeqIO.parse(open(filename), "fasta")

    seq_dict = dict()
    for record in seqs:
        seq_dict[record.id] = record

    return seq_dict

def read_seq_table(seq_table, ref_positions):
    """
    Parameters
    ----------
    seq_table : a Pandas data frame
        a variation table created by consensus3.py
    ref_positions : a list
        a list of columns in seq_table to be checked for nucleotides

    Returns
    -------
    seq_ids : a list
        a list of sequence IDs containing valid nucletide values at ref_positions 

    """
    df = pd.read_csv(seq_table, header=0)
    
    cols = ['Pos ' + str(x) for x in ref_positions]
    
    seq_ids = []
    for _, row in df.iterrows():
        ok = True
        for col in cols:
            if not row[col] in ['A', 'C', 'G', 'T']:
                ok = False
                break
        if ok:
            seq_ids.append(row['ID'])
            
    return seq_ids

def remove_duplicates(seqs, id_list):
    seq_dict = dict()
    for id in id_list:
        if not seqs[id].seq in seq_dict:
            seq_dict[seqs[id].seq] = seqs[id]

    return seq_dict

def main():
    date = '2020_10_07'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    in_file = base + 'sequences_no_dups.fasta'
    out_file = base + 'sequences_filtered.fasta'
    seq_table = base + 'sars_cov_2_variation_ncbi_no_dups_98.0.csv'
    ref_positions = [241, 3037, 14408, 23403]
    
    os.chdir(base)
        
    id_list = read_seq_table(seq_table, ref_positions)
    all_seqs = read_fasta_file(in_file)

    seqs = remove_duplicates(all_seqs, id_list)

    with open(out_file, 'w') as f:
        for s, seq_rec in seqs.items():
            SeqIO.write(seq_rec, f, 'fasta')


if __name__ == "__main__":
    main()
