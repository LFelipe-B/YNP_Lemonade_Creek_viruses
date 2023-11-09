#!/usr/bin/python


# Usage: python Get_viral_scaffold_approx_class_V3_publication.py *_blastp_RefSeqComplete.outfmt6.gz_TOPHITS_edit.tsv

# codes mostly from here:
# codes from here: https://stackoverflow.com/questions/33290374/pandas-pivot-table-column-names and here: https://sparkbyexamples.com/pandas/pandas-remove-duplicate-columns-from-dataframe/
# diamond blastp ouput tabled edited to include in column 1 the scaffold identification

# import lib
import pandas as pd
import sys
import csv

# open input file
f = open(sys.argv[1])

# read input file
data = pd.read_csv((f), delimiter='\t',  index_col = 0, keep_default_na=False) # or data = pd.read_csv('CreekBiofilm_2_diamond_vir.txt_top_hits_headers.tsv', delimiter='\t', index_col = 0, keep_default_na=False)

### First part of the code -> get viral counts and print True/False if is a viral contig

# String to be searched
search ="Viruses"

# count occurrence distribution of string (word) and create new column (count_tax, search=virus) based on species name
data["count_virus"]= data["sskingdoms"].str.count(search)

# create pivot dataframe (df) in this case "pivot" with 3 cols scaffold_name, count (number of lines (hits) of Query), and sum (number of times of searched word - string was found in this case "Viruses")
pivot = data.pivot_table(index= 'scaffold_name', values= "count_virus", aggfunc={'count','sum'})

# create new col in pivot df with name "calc" and get results from number of lines of number of virus hits divided by number of hits in percentage
pivot['calc'] = pivot['sum']/pivot['count']*100

# create new col in pivot df with name "Viral_contig" that print" true" if equal or greather than 60% or" false" if less
pivot['Viral_contig'] = pivot['calc'] >= 60

# change col names from 'count, sum, next, keep/discard' to 'n of hits', 'n of viruses','half n of hits', 'Taxon dist HGT'. This should be raw table#
pivot.columns = ['N_hits', 'N_viruses','Percent_viruses_60', 'Viral_contig']

# reset the index and save as the same object
pivot = pivot.reset_index()

# print result of final dataframe removing the index (integer numeric column)
pivot.to_csv(sys.argv[1] + '_viral_scaffold_approx_class.csv', index=False)