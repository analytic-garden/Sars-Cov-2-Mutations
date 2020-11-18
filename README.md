See https://analyticgarden.blogspot.com/2020/11/sars-cov-2-co-varying-changes.html
for examples.

Pipeline

Python

1) Download SARS-CoV-2 FASTA and GenBank sequences from
https://www.ncbi.nlm.nih.gov/datasets/coronavirus/genomes/
unzip into ncbi_dataset
data is in ncbi_dataset/data/genomic.fna and ncbi_dataset/data/genomic.gbff

2) Remove sequences with bad dates
run remove_dups_dates.py
Produces new FASTA file sequences_valid_dates.fasta and NC_045512.fasta

3) Align sequences
Align sequences to reference sequence NC_045512.2
time ~/anaconda3/bin/mafft --auto --preservecase --addfragments sequences_valid_dates.fasta NC_045512.fasta > sequences_valid_dates_aln.fasta
or
time ~/anaconda3/bin/mafft --auto --preservecase --thread -1 --addfragments sequences_valid_dates.fasta NC_045512.fasta > sequences_valid_dates_aln.fasta

Produces sequences_valid_dates_aln.fasta

4)
edit consensus cutoff and min_qual_cutoff in each of these to change
the number of columns used from the alignment.

run variation.py 
produces: sars_cov_2_variation_ncbi_valid_dates_98.0.csv
positions with signifcant variation by NCBI ID
98.0 is the consensus cutoff

run mutations.py
produces: sars_cov_2_ncbi_ncbi_mut_valid_dates.csv
positions with significant mutations

run MI.py
produces: MI_ncbi_valid_dates.csv
mutual information among aligned positions

R

5) Plot MI graph
mi_2020_11_11 <- read.csv('MI_ncbi_valid_dates.csv')
plot_MI(mi_2020_11_11)
2020_11_11 is the download date. The object name can be anything.

6) Read sars_cov_2_variation_ncbi_valid_dates_98.0.csv into dataframe
cv_2020_11_11 <- read.csv('sars_cov_2_variation_ncbi_valid_dates_98.0.csv')
# set Collection.Date to datetime
cv_2020_11_11$Collection.Date <- as.Date(cv_2020_11_11$Collection.Date)

7) run plot_varying_pct3(cv_2020_11_11, positions = c(241, 3037, 14408, 23403))
to get plots for groups of sites
run plot_varying_nucs2(cv_2020_11_11,positions = c(241, 3037, 14408, 23403))
to get bar plots of counts
c(241, 3037, 14408, 23403) is a list of position you wish to plot

11) Plot combos
cv_table_241 <- covary_table2(cv_2020_11_11, positions = c(241, 3037, 14408, 23403))
ggplot(cv_table_241, aes(x = Collection.Date, y=Count, fill=Nucleotides)) +
		     geom_bar(stat='identity', width=1)
		     
get counts of each combination
temp <- cv_table_241 %>%
     	   group_by(Nucleotides) %>%
     	   select(Count, Nucleotides) %>%
	   summarise_at(c("Count"), sum)

