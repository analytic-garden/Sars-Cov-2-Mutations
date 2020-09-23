"""
reformat_aln_headers.py - add date and location to FASTA header
author: Bill Thompson
license: GPL 3
copyright: Apr. 8, 2020
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from consensus3 import read_alignment_file, read_genbank_file, format_date

def get_quartet(df, id, seq_positions):
    """
    get_quartet - get the 4 mnuclotide positions from GSAID clade GH

    Parameters
    ----------
    df : str
        a data frame of the variation file created by position3.py
    id : str
        a GenBank sequence ID
    seq_positions : list
        list of variant positions to retreive

    Returns
    -------
    seq_str : str
        the formatted nucleotides at seq_position

    """
    seq_str = '-'.join(df.loc[id, seq_positions].to_list())

    return seq_str

def main():
    date = '2020_08_28'
    path = '/mnt/g/Covid-19/' + date + '/'
   
    align_file = path + 'sequences_filtered_aln_red.fasta'
    genbank_file = path + 'sequences.gb'
    out_file = path + 'sequences_filtered_aln2_red.fasta'
    csv_file = path + 'sars_cov_2_variation_ncbi_no_dups_98.0.csv'
    
    seq_positions = ['Pos 241', 'Pos 3037', 'Pos 14408', 'Pos 23403']

    df = pd.read_csv(csv_file, header=0, index_col=0)

    align = read_alignment_file(align_file)
    gb = read_genbank_file(genbank_file)

    with open(out_file, 'w') as f:
        for id, seq in align.items():
            gb_seq = gb[id]
            if 'country' not in gb_seq.features[0].qualifiers:
                country = 'Unknown'
            else:
                country = gb_seq.features[0].qualifiers['country'][0]
            date = gb_seq.features[0].qualifiers['collection_date'][0]
            dt = format_date(date)
            var1 = get_quartet(df, id, seq_positions)
            seq_rec = SeqRecord(seq.seq,
                                id = id,
                                description = ' '.join([country,
                                                         dt.strftime('%Y-%m-%d'),
                                                         var1]))
            SeqIO.write(seq_rec, f, 'fasta')


if __name__ == "__main__":
    main()
