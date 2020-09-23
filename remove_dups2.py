"""
remove_dupss.py - remove duplicate sequences and 
                  add date and location to FASTA header
author: Bill Thompson
license: GPL 3
copyright: Apr. 30, 2020
"""

import os
from Bio import SeqIO

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


def remove_duplicates(seqs, ref_id):
    """
    remove_duplicates - remove duplicate sequences
    
    Running this will give a warning. We are comparing the sequences as 
    strings.

    Parameters
    ----------
    seqs : dictionary
        BioPython SeqRecs read by read_fasta_file
        key is sequence ID
    ref_id : str
        ID string of reference sequence

    Returns
    -------
    seq_dict : dictionary
        a dictionary with sequence as key and SeqRec as value

    """
    seq_dict = dict()

    seq_dict[seqs[ref_id].seq] = seqs[ref_id]

    for id, seq_rec in seqs.items():
        if not seq_rec.seq in seq_dict:
            seq_dict[seq_rec.seq] = seq_rec

    return seq_dict


def main():
    path = '/mnt/g/Covid-19/2020_09_04/'
    
    fasta_file = 'sequences.fasta'
    out_file = 'sequences_no_dups.fasta'
    ref_id = 'NC_045512.2'
    
    os.chdir(path)
    
    all_seqs = read_fasta_file(fasta_file)

    seqs = remove_duplicates(all_seqs, ref_id)

    with open(out_file, 'w') as f:
        for s, seq_rec in seqs.items():
            SeqIO.write(seq_rec, f, 'fasta')


if __name__ == "__main__":
    main()
