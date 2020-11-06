"""
consensus3.py - format and save SARS-CoV-2 data from NCBI
author: Bill Thompson
license: GPL 3
copyright: Mar. 27, 2020
"""

import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from collections import defaultdict
from sklearn import metrics
import numpy as np
from datetime import datetime
import sys

"""
MI - calculate the mutual information between two columns of a 
     multisequence alignment

arguments:
col1, col2 - two columns from the alignment
pseudo - pseudocounts

returns:
the mutual information in bits

requires:
len(col1) == len(col2)
pseudo >= 0
np.log2
"""

def MI(col1, col2, pseudo = 0.05):
    return metrics.mutual_info_score(col1, col2)/np.log(2)

"""
MI_table - generate a table of mutual information values from all paires of
           columns in a  dictionary of aligned nucleotides

arguments:
variant_cols - a dictionary of columns. 
               key = aligned column number, value = column data
pos_map - a dictionary, key = aligned column position, value = reference coulum  position
pseudo - pseudocounts

returns:
a dictionary -  sorted by value
     key = string "col1_number, col2_number", 
     value = mutual information between columns

requires:
pseudo >= 0
"""

def MI_table(variant_cols, pos_map, pseudo = 0.05):
    var_cols = list(variant_cols.keys())
    mi_tab = dict()
    for i in range(len(var_cols)-1):
        for j in range(i+1, len(var_cols)):
            k = str(pos_map[var_cols[i]]+1) + ', ' + str(pos_map[var_cols[j]]+1)
            mi_tab[k] = MI(variant_cols[var_cols[i]][1],
                           variant_cols[var_cols[j]][1],
                           pseudo = pseudo)
    mi_tab = {k: v for k, v in sorted(mi_tab.items(), key=lambda item: item[1], reverse=True)}

    return mi_tab


"""
read_alignment_file - read a collection of aligned sequences in FASTA format

arguments:
filename - the name of the FASTA file

returns:
a dictionary. The keys are sequence IDs. The values are BioPython SeqRecords.
"""
def read_alignment_file(filename):
    alignment = AlignIO.read(open(filename), "fasta")

    align = dict()
    for record in alignment:
        align[record.id] = record

    return align


"""
read_genbank_file - read a collection of GenBank records downloaded from NCBI

arguments:
filename - the name of the GenBank file

returns:
a dictionary. The keys are sequence IDs. The values are BioPython SeqRecords.
"""
def read_genbank_file(filename):
    records = SeqIO.parse(open(filename), "genbank")

    gb = dict()
    for record in records:
        gb[record.id] = record

    return gb

"""
ref_pos_to_alignment - create a map of aligned positions to reference 
                       sequence positions

arguments:
align_file - FASTA file of aligned sequences
ref_id - ID of the reference sequence

returns:
a dictionary, key = aligned column position, value = reference coulum  position
"""

def ref_pos_to_alignment(align_file, ref_id):
    algn = read_alignment_file(align_file) 
    ref_seq = algn[ref_id].seq

    pos_map = dict()
    ref_pos = 0    # 0 based
    for pos in range(len(ref_seq)):
        if ref_seq[pos] in ['A', 'C', 'G', 'T']:
            pos_map[pos] = ref_pos
            ref_pos += 1
        else:
            pos_map[pos] = -1

    return pos_map
                    
"""
get_varying_columns - find columns in alignment that vary by more than a
                      given amount

arguments:
align - a Bio.Align object produced by Bio.AlignIO
consensus_cutoff - upper value for variation in an alignment column
start, end - starting and ending columns in the alignment

returns:
a dictionary - key = a column number. Thse colummns have identity <= cutoff
               value = tuple (percent of variation, 
                              a list containing the aligned column nucleotides)

requires:
consensus_cutoff, start, end >= 0
start, end with alignment column bounds
"""
def get_varying_columns(align,
                        start = 130, end = 29840):
    variant_cols = dict()
    for col in range(start, end):
        # c = Counter(align[:, col]).most_common(1)
        c = Counter(align[:, col])
        common = c.most_common(1)
        if common[0][0] in ['A', 'C', 'G', 'T']:
            denom = sum([c[k] for k in ['A', 'C', 'G', 'T', '-']])
            pct = common[0][1] / denom
            if pct < 1:
                variant_cols[col] = (pct, list(align[:, col]))

    return variant_cols

"""
format_date - a helper function. Dates from GenBank file are not uniformly formatted

arguments:
date - a string containing a date

returns:
a datetime date
"""
def format_date(date):
    # try and guess the format
    try:
        dt = datetime.strptime(date, '%d-%b-%Y')
    except ValueError:
        try:
            dt = datetime.strptime(date, '%b-%Y')
        except ValueError:
            try:
                dt = datetime.strptime(date, '%m/%d/%Y')
            except ValueError:
                try:
                    dt = datetime.strptime(date, '%Y-%m-%d')
                except:
                    try:
                        dt = datetime.strptime(date, '%Y-%m')
                    except:
                        dt = datetime.strptime(date, '%Y')

    return dt

"""
init_dataframe - a helper funtion to initialize the first few columns of
                 csv_file_base 

arguments:
genbank - a dictionary, key = GenBank ID, value - Bio.SeIO GenBank record
align - a collection of alignment records read by Bio.AlignIO

returns:
a Pandas dataframe with columns ID, Collection Date, Country of Origin
"""

def init_dataframe(genbank, align):
    # store the data in a Pandas dataframe
    # BioPython uses 0-based indexing. GenBank uses 1-based indexing
    df = pd.DataFrame(columns = ['ID', 'Collection Date', 'Country'])
    for seq_rec in align:
        id = seq_rec.id
    
        # get country and date info.
        gb = genbank[id]
        if 'country' not in gb.features[0].qualifiers:
            country = 'Unknown'
        else:
            country = gb.features[0].qualifiers['country'][0]
        date = format_date(gb.features[0].qualifiers['collection_date'][0])

        df = df.append({'ID': id,
                        'Collection Date': date,
                        'Country': country},
                        ignore_index=True)

    return df

"""
format_qual - format GenBank records qulaifiers.
              GenBank qualifiers contain commas. Since we are writing to a CVS
              we surround strings containing commas with quotes.
              See BioPython SeqIO docs for info about qualifiers.

arguments:
qualifiers - a dictionary, key = GenBank ID, value =  GenBank qualifier
id - a GenBank ID

returns:
a string - the qualifier string, possibly surrounded by quotes, or an empty string
"""
def format_qual(qualifiers, id):
    if not id in qualifiers:
        return ''

    qual = qualifiers[id][0]
    if ',' in qual:
        qual = '"' + qual + '"'

    return qual
                         
"""
count_mutations - count the variation in an alignment
                  This is a big ugly function that counts variation in columns
                  with more variation than consensus_cutoff. It returns a
                  collection  of data for the mutation CSV file.

arguments:
align - a Bio.Align object produced by Bio.AlignIO
pos_map - a dictionary, key = aligned column position, value = reference coulum  position
ref_seq - the reference sequence. A Bio.Seq object
consensus_cutoff - upper value for variation in an alignment column
start, end - starting and ending columns in the alignment

returns:
a dictionary
alignment columns, reference columns - list of column numbers in alignment and reference sequence
A,C,G,T,insertions - list the count of nucleutides and insertion in aligned columns
other - list of the counts of everything else (N's ambiguous calls, etc.)
consensus - a list consensus nucleotide for each column
consensus %' - a list consensus percent of consensus nucleotide for each column
alt - a list of the next popular nucleotide in each column
codons - list of codons containing consensus nucleotide
aas - list amino acid of codon containing consensus nucleotide
alt_codons - a list of the codons containing the enxt popular nucleotide 
alt_aas - a list of the amino acids of the codons containing the enxt popular nucleotide 
products, notes - protein and info GenBank listed by GenBank
"""

def count_mutations(align, pos_map, ref_seq,
                    start=130, end=29840, consensus_cutoff=0.8):
    rows = align.__len__()

    align_cols = []
    ref_cols = []
    nts = {'A': [], 'C': [], 'G': [], 'T': [], '-': []}
    other = []
    consensus = []
    consensus_pct = []
    alt_consensus = []
    
    # run through columns and find coulmns with significant variation
    for col in range(start, end):          
        c = Counter(align[:, col])
        if len(c) > 1:
            common = c.most_common(2)
            
            if common[0][0] in ['N', '-']:
                continue   # ignore colums that are mostly N's or gaps
            
            if len(c) >= 2:
                if common[1][0] in ['N', '-']:  # ignore N's
                    continue

            # count nucleotides
            counts = defaultdict(int)
            for nt, count in c.items():
                counts[nt] += count

            for nt in nts.keys():
                nts[nt].append(counts[nt])
            other.append(rows - sum([counts[k] for k in nts.keys()]))

            # get the consensu and next best
            consens = common[0][0]
            consensus.append(consens)
            try:
                consensus_pct.append((counts[consens] / sum([counts[k] for k in nts.keys()])) * 100)
            except:
                print('Invalid column:', col)
                sys.exit(1)                
            alt_consensus.append(common[1][0])
            
            align_cols.append(col+1)
            ref_cols.append(pos_map[col]+1)

    codons = []
    aas = []
    alt_codons= []
    alt_aas = []
    notes = []
    products = []
    refs = []
    aa_pos = []
    feat_type = []
    mutations = []
    new_nucs = []
    # find the codons and amino acids
    for position, pct, alt, consens in zip(ref_cols,
                                           consensus_pct,
                                           alt_consensus,
                                           consensus):
        if pct/100 < consensus_cutoff: # only use significan columns
            seq_info = find_location_feature(ref_seq, position-1, alt, consens)
            aas.append(str(seq_info['aa']))
            feat_type.append(seq_info['feature_type'])
            
            if seq_info['feature_type'].find('UTR') == -1:
                codons.append(str(seq_info['codon']))
                alt_codons.append(str(seq_info['alt_codon']))
                alt_aas.append(str(seq_info['alt_aa']))
                mutations.append(str(seq_info['aa']) + \
                                 str(int(seq_info['aa_pos'])+1) + \
                                 str(seq_info['alt_aa']))
                new_nucs.append(str(seq_info['new_nuc']))
            else:
                codons.append('')
                alt_codons.append('')
                alt_aas.append('')
                mutations.append('')
                new_nucs.append('')
                
            notes.append(format_qual(seq_info['qualifiers'], 'note'))
            products.append(format_qual(seq_info['qualifiers'], 'product'))
            refs.append(seq_info['ref_nuc'])
            aa_pos.append(seq_info['aa_pos']+1)
        else:
            codons.append('')
            aas.append('')
            alt_codons.append('')
            alt_aas.append('')
            notes.append('')
            products.append('')
            refs.append('')
            aa_pos.append('')
            feat_type.append('')
            mutations.append('')
            new_nucs.append('')

    return {'alignment columns': align_cols,
            'reference positions': ref_cols,
            'A': nts['A'],
            'C': nts['C'],
            'G': nts['G'],
            'T': nts['T'],
            'insertions': nts['-'],
            'other': other,
            'ref': refs,
            'consensus': consensus,
            'consensus %': consensus_pct,
            'alt': new_nucs,
            'feature_type': feat_type,
            'codons': codons,
            'ref_aa': aas,
            'alt_codons': alt_codons,
            'alt_aas': alt_aas,
            'aa_pos': aa_pos,
            'mutations': mutations,
            'products': products,
            'notes': notes
            }

"""
find_location_feature - get features at position from reference sequence

arguments:
ref_seq - a BioSeq GenBank record
position - position in sequence
alt - an alternate nucleotide, most likely the next best nucleotide of consensus

returns:
a dictionary
conserved_position - the input position,
start, end - start and end of feature
feature_trype - a string, the feature type gene, UTR etc.
codon_start, codon_end - atrt and end of codon
codon - the codon,
aa - the aa at position
alt_codon, alt_aa - codon and amino acid if alt is substituted at position
seq - sequence of feature
qualifiers - feature qualifiers See Bio.Seq docs for details

requires:
position - must be in range of reference sequence length
alt - A, C, G, or T

"""

def find_location_feature(ref_seq, position, alt, consens):
    feat = ''
    for f in ref_seq.features:
        if position in f:
            feat = f
            
    feat_type = feat.type
       
    for p in feat.location.parts:
        start = p.start.real
        end = p.end.real
        if position >= start and position < end:
            break;
             
    sequence = ref_seq.seq[start:end]
        
    # start = feat.location.start.real
    # end = feat.location.end.real
    # sequence = ref_seq[start:(end+1)]
                       
    offset = position - start
    codon_pos = offset % 3
    codon_start = offset - codon_pos
    codon_end = codon_start + 2
    # codon = feat.extract(ref_seq).seq[codon_start:(codon_end+1)]
    codon = sequence[codon_start:(codon_end+1)]
    aa = codon.translate()
    
    ref_nuc = sequence[offset]
    new_nuc = alt
    if alt == ref_nuc:
        new_nuc = consens
    aa_pos = codon_start / 3
    
    seq_list = list(str(sequence))
    seq_list[offset] = new_nuc
    new_seq = SeqRecord(Seq(''.join(seq_list)))
    alt_codon = new_seq.seq[codon_start:(codon_end+1)]
    alt_aa = alt_codon.translate() if new_nuc in ['A', 'C', 'G', 'T'] else ''
    
    return {'conserved_position': position,
            'start': start,
            'end': end,
            'feature_type': feat_type,
            'codon_start': codon_start,
            'codon_end': codon_end,
            'codon': codon,
            'aa': aa,
            'alt_codon': alt_codon,
            'alt_aa': alt_aa,
            'ref_nuc': ref_nuc,
            'aa_pos': aa_pos,
            'seq': sequence,
            'qualifiers': feat.qualifiers,
            'new_nuc': new_nuc}

"""
count_variation_from_ref - count the differences from the reference sequence

arguments:
align - a Bio.Align object produced by Bio.AlignIO
align_file - file containing alignmet in FASTA format
ref_id - the ID of the refeence sequence
start, end - start and end positions in alignmet to compare

returns:
a list of the differences from the reference sequence for each 
sequences in the alignment

"""

def count_variation_from_ref(align, align_file, ref_id,
                             start = 130, end = 29840):
    align_dict = read_alignment_file(align_file)
    ref_seq = align_dict[ref_id]

    diffs = []
    for seq_rec in align:
        if seq_rec.id == ref_id:
            print()
        seq = seq_rec.seq
        diff_count = sum(seq[i] != ref_seq.seq[i] for i in range(start, end))
        diffs.append(diff_count)

    return diffs

def get_col_range(align, min_quality = 0.99):
    """

    Parameters
    ----------
    align : a BioPython Bio.Align.MultipleSeqAlignment object
    min_quality : a float
        the minimum column quality
        quality is the percent of A, C, G, T in column

    Returns
    -------
    start, end : ints
        the start and ending positions of positions with sufficient quality

    """    
    start = 0
    for col in range(0, align.get_alignment_length()):
        c = Counter(align[:, col]).most_common(4)
        qual = sum([x[1] for x in c if x[0] in ['A', 'C', 'G', 'T']])/len(align[:, col])
        if qual > min_quality:
            start = col
            break
        
    end = 0
    for col in range(align.get_alignment_length()-1, -1, -1):
        c = Counter(align[:, col]).most_common(4)
        qual = sum([x[1] for x in c if x[0] in ['A', 'C', 'G', 'T']])/len(align[:, col])
        if qual > min_quality:
            end = col
            break
    
    return start, end
    
def main():
    date = '2020_10_07'
    base = '/mnt/g/Covid-19/' + date + '/' 
    
    align_file =  base + 'sequences_no_dups_aln.fasta'
    genbank_file = base + 'seqs.gb/ncbi_dataset/data/genomic.gbff'
    csv_file_base = base + 'sars_cov_2_variation_ncbi_no_dups_'
    mutation_csv_base = base + 'sars_cov_2_ncbi_ncbi_mut_no_dups_'
    mutual_info_csv = base + 'MI_ncbi_no_dups.csv'

    ref_id = 'NC_045512.2'

    # change these for different alignments and mutation level cutoffs
    min_col_quality = 0.90
    consensus_cutoff = 0.98
    mi_cutoff = 0.5
    csv_file = csv_file_base + str(consensus_cutoff * 100) + '.csv'
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

    df = init_dataframe(genbank, align)
        
    # get varying columns
    variant_cols = get_varying_columns(align, 
                                       consensus_cutoff = consensus_cutoff,
                                       start = start, end = end)

    diffs = count_variation_from_ref(align, align_file, ref_id, start, end)
    df['Diffs from Ref'] = diffs
    
    # wite the data
    print('Variant Positions')
    for col in variant_cols:
        df['Pos ' + str(pos_map[col]+1)] = variant_cols[col][1]
        print('Pos ' + str(pos_map[col]+1) + ':', variant_cols[col][0])
    print()

    df = df.sort_values('Collection Date')
    df.to_csv(csv_file, index=False)

    # mutual information
    print('Mutual Infomation')
    mi_tab = MI_table(variant_cols, pos_map)
    mi_dict = {'Position_1': [], 'Position_2': [], 'MI': []}
    for k,mi in mi_tab.items():
        p1, p2 = k.split(',')
        mi_dict['Position_1'].append(int(p1))
        mi_dict['Position_2'].append(int(p2))
        mi_dict['MI'].append(mi)
        if mi >= mi_cutoff:
            print('MI(', k,  ') =', mi)
    df_mi = pd.DataFrame(mi_dict)
    df_mi.to_csv(mutual_info_csv, index=False)

    # mutations = count_mutations(align, pos_map, ref_seq, 
    #                             start = start, end = end,
    #                             consensus_cutoff=consensus_cutoff)
    # df2 = pd.DataFrame(mutations)
    # df2 = df2.sort_values('consensus %')
    # df2.to_csv(mutation_csv_file, index=False)

if __name__ == "__main__":
    main()
