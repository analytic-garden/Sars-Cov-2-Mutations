"""
remove_dups_dates.py - remove duplicate sequences and sequences with
                       invalid collection dates.
author: Bill Thompson
license: GPL 3
copyright: Apr. 30, 2020
"""

from Bio import SeqIO
from datetime import datetime
from consensus3 import read_genbank_file

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
        second term is TRUE if date is complte
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

def remove_duplicates(seqs, gb, ref_id):
    """
    remove_duplicates - remove duplicate and sequences with invalid 
                        dates 
    
    Running this will give a warning. We are comparing the sequences as 
    strings.

    Parameters
    ----------
    seqs : dictionary
        BioPython SeqRecs read by read_fasta_file
        key is sequence ID
    gb : dictionary
        BioPython SeqRecs read by read_genbank_file
        key is sequence ID
    ref_id : str
        ID string of reference sequence

    Returns
    -------
    seq_dict : dictionary
        a dictionary with sequence as key and SeqRec as value

    """
    seq_dict = dict()
    dup_seqs = list()
    bad_dates = list()

    seq_dict[seqs[ref_id].seq] = seqs[ref_id]

    for id, seq_rec in seqs.items():
        gb_rec = gb[id]
        collection_date = gb_rec.features[0].qualifiers['collection_date'][0]
        formatted_date, ok = check_date(collection_date)

        if ok: 
            if not seq_rec.seq in seq_dict:
                seq_dict[seq_rec.seq] = seq_rec
            elif id != ref_id:
                dup_seqs.append(seq_rec.id)
        else:
            bad_dates.append(seq_rec.id)

    return seq_dict, dup_seqs, bad_dates

def main():
    path_date = '2020_10_27'
    base = '/mnt/g/Covid-19/' + path_date + '/' 

    fasta_file = base + 'ncbi_dataset/data/genomic.fna'
    genbank_file = base + 'ncbi_dataset/data/genomic.gbff'
    refseq_file = base + 'NC_045512.fasta'
    out_file = base + 'sequences_no_dups.fasta'
    ref_id = 'NC_045512.2'
        
    gb = read_genbank_file(genbank_file)
   
    all_seqs = read_fasta_file(fasta_file)

    seqs, dup_seqs, bad_dates = remove_duplicates(all_seqs, gb, ref_id)
    
    f1 = open(out_file, 'w')
    f2 = open(refseq_file, 'w')
    for s, seq_rec in seqs.items():
        if seq_rec.id == ref_id:
            SeqIO.write(seq_rec, f2, 'fasta')
        else:
            SeqIO.write(seq_rec, f1, 'fasta')
    f1.close()
    f2.close()

if __name__ == "__main__":
    main()
