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


def main():
    id_file = 'sequences.acc'
    out_file = 'sequences.gb'
    
    os.chdir('/mnt/g/Covid-19/2020_09_04')

    id_list = read_id_file(id_file)

    Entrez.email = 'analytic_garden@gmail.com'
    
    out_handle = open(out_file, 'w')
    start = 0
    ret_max = 500
    while start < len(id_list):
        end = min(start+ret_max, len(id_list))
        in_handle = Entrez.efetch(db="nucleotide", 
                                  id=id_list[start:end], 
                                  rettype="gb")
        recs = in_handle.read()
        out_handle.write(recs)
        print('records', start, 'to', end-1)
        start += (end - start)
    out_handle.close()
    in_handle.close()
    
if __name__ == '__main__':
    main()
