"""
fetch_genbank.py - retreive sequences from GenBank
author: Bill Thompson
license: GPL
copyright: Mar. 15, 2020
"""

import os
from Bio import Entrez

def read_id_file(filename):
    """
    read_id_file - read a file of accession numbers

    Parameters
    ----------
    filename : str
        name of the file, sequences.acc.

    Returns
    -------
    ids : list
        list of accession numbers

    """
    with open(filename) as f:
        ids = f.readlines()

    ids = [x.strip() for x in ids]

    return ids

def find_missing_recs(id_list, recs):
    """
    find_missing_recs - get a list of GenBank records not downloades

    Parameters
    ----------
    id_list : list
        a list of IDs intended for download from GenBank.
    recs : string
        a large string containing GenBank records.

    Returns
    -------
    missing : list
        a list of GenBank IDs not returned.

    """
    missing = []
    for id in id_list:
        if recs.find(id) == -1:
            print('Missing GenBank record:', id)
            missing.append(id)
            
    return missing

def main():
#    id_file = 'sequences.acc'
    id_file = 's.acc'
    out_file = 'sequences.gb'
    
    os.chdir('/mnt/g/Covid-19/2020_10_07')

    id_list = read_id_file(id_file)

    Entrez.email = 'analytic_garden@gmail.com'
    
    missing_list = list()
    
    out_handle = open(out_file, 'w')
    start = 0
    ret_max = 500
    while start < len(id_list):
        end = min(start+ret_max, len(id_list))
        in_handle = Entrez.efetch(db="nucleotide", 
                                  id=id_list[start:end], 
                                  rettype="gb")
        recs = in_handle.read()
        if recs.count('LOCUS') != len(id_list[start:end]):
            missing = find_missing_recs(id_list[start:end], recs)
            missing_list += missing
        
        out_handle.write(recs)
        print('records', start, 'to', end-1, 'count:', recs.count('LOCUS'))
        start += (end - start)
    out_handle.close()
    in_handle.close()
    
    print('Missing IDs')
    for id in missing_list:
        print(id)
    
if __name__ == '__main__':
    main()
