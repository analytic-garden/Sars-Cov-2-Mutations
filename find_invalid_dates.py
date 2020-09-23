#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
find_invalid_dates.py

@author: Bill Thompson
@license: GPL 3
@copyright: Sept. 14, 2020
"""

from datetime import datetime
from consensus3 import read_genbank_file
from remove_dups_dates import read_fasta_file

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

def main():
    date = '2020_09_04'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    fasta_file = base + 'sequences_no_dups.fasta'
    genbank_file = base + 'sequences.gb'
    
    all_seqs = read_fasta_file(fasta_file)
    gb = read_genbank_file(genbank_file)
    
    for id, record in gb.items():
        collection_date = record.features[0].qualifiers['collection_date'][0]
        formatted_date, ok = check_date(collection_date)
        if (not ok) and id in all_seqs :
            print(id, collection_date, formatted_date)

if __name__ == "__main__":
    main()

