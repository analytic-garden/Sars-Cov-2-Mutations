"""
remove_bad_dates.py - sequences with invalid collection dates.
    produces two files: 
        NC_045512.fasta containing the reference sequence
        sequences_valid_dates.fasta containging all thes equences with valid dates        
author: Bill Thompson
license: GPL 3
copyright: Oct. 29, 2020
"""

from Bio import SeqIO
from datetime import datetime
from consensus3 import read_genbank_file

def read_fasta_file(filename):
    """
    read_fasta_file - read a collection sequences in FASTA format

    arguments:
        filename - the name of the FASTA file

    returns:
        a dictionary. The keys are sequence IDs. The values are BioPython SeqRecords.
    """
    seqs = SeqIO.parse(open(filename), "fasta")

    seq_dict = dict()
    for record in seqs:
        seq_dict[record.id] = record

    return seq_dict

def check_date(date):
    """
    check_date - check for complete date, month, day, year

    Parameters
    ----------
    date : str
        a GenBank collection date.

    Returns
    -------
    dt : tuple
        (datetime object, boolean) 
        second term is TRUE if date is complete
    """
    # try and guess the format
    date_ok = True
    try:
        dt = datetime.strptime(date, '%d-%b-%Y')
    except ValueError:
        try:
            dt = datetime.strptime(date, '%b-%Y')
            date_ok = False
        except ValueError:
            try:
                dt = datetime.strptime(date, '%m/%d/%Y')
            except ValueError:
                try:
                    dt = datetime.strptime(date, '%Y-%m-%d')
                except:
                    try:
                        dt = datetime.strptime(date, '%Y-%m')
                        date_ok = False
                    except:
                        dt = datetime.strptime(date, '%Y')
                        date_ok = False

    return (dt, date_ok)

def main():
    path_date = '2020_11_11'
    base = '/mnt/g/Covid-19/' + path_date + '/' 

    fasta_file = base + 'ncbi_dataset/data/genomic.fna'
    genbank_file = base + 'ncbi_dataset/data/genomic.gbff'
    refseq_file = base + 'NC_045512.fasta'
    out_file = base + 'sequences_valid_dates.fasta'
    ref_id = 'NC_045512.2'
        
    gb = read_genbank_file(genbank_file)   
    all_seqs = read_fasta_file(fasta_file)

    with open(refseq_file, 'w') as f:
        SeqIO.write(all_seqs[ref_id], f, 'fasta')
        
    seqs = dict()        
    for seq_id, seq_rec in all_seqs.items():
        gb_rec = gb[seq_id]
        collection_date = gb_rec.features[0].qualifiers['collection_date'][0]
        formatted_date, ok = check_date(collection_date)
        if ok:
            seqs[seq_id] = seq_rec

    with open(out_file, 'w') as f1:
        for seq_id, seq_rec in seqs.items():
            SeqIO.write(seq_rec, f1, 'fasta')

if __name__ == "__main__":
    main()
