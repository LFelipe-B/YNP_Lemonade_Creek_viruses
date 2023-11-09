# usage
python get_most_freq_tax_consensus_publication.py Creek_1_all_vir_scaffold_final_classification_edited_example_file

# upload python library
import pandas as pd
import sys
import csv

# import and load table file
f = open(sys.argv[1])

# upload input table (use two columns only "Metag_scaffold" and "Lineage")
data = pd.read_csv((f), delimiter='\t', index_col = 0)

# group by Bin, and print most common value https://stackoverflow.com/questions/15222754/groupby-pandas-dataframe-and-select-most-common-value
data_grouped = data.groupby(['Metag_scaffold'], sort=True)['Lineage'].agg(pd.Series.mode)

# print final table
data_grouped.to_csv(sys.argv[1] + '_consensus.tsv')
